#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// get car index in a given lane which is ahead of a given s
int getFrontCar(int lane,  double s, vector<vector<double>> sensor_fusion) {

  int idx = -1;
  double closest_dis = 6945.554;
  //find if there is a car ahead and it is too close
  for(int i=0; i < sensor_fusion.size(); i++) {

    //car in current lane
    double d = sensor_fusion[i][6];
    if (d >= lane*4. && d <= (lane+1)*4.) {
      double other_car_s = sensor_fusion[i][5];
      //if other car is ahead of ego and closer then previous found car
      if (other_car_s > s && (other_car_s-s) < closest_dis) {
        closest_dis = other_car_s-s;
        idx = i;
      }
    }
  }
  return idx;
}

// get car index in a given lane which is behind  of a given s
int getBackCar(int lane,  double s, vector<vector<double>> sensor_fusion) {

  int idx = -1;
  double closest_dis = 6945.554;
  //find if there is a car ahead and it is too close
  for(int i=0; i < sensor_fusion.size(); i++) {

    //car in current lane
    double d = sensor_fusion[i][6];
    if (d >= lane*4. && d <= (lane+1)*4.) {
      double other_car_s = sensor_fusion[i][5];
      //if other car is ahead of ego and closer then previous found car
      if (other_car_s < s && (s-other_car_s) < closest_dis) {
        closest_dis = s-other_car_s;
        idx = i;
      }
    }
  }
  return idx;
}

bool safetoChangeTo(int new_lane, int lane, int prev_size, double car_s, double speed, vector<vector<double>> sensor_fusion) {
  int f_idx = getFrontCar(new_lane, car_s, sensor_fusion);
  int b_idx = getBackCar(new_lane, car_s, sensor_fusion);
  bool front_clear = true;

  if (f_idx != -1) {
    double vx = sensor_fusion[f_idx][3];
    double vy = sensor_fusion[f_idx][4];
    double other_car_speed = sqrt(vx*vx + vy*vy);
    double other_car_s = sensor_fusion[f_idx][5];

    //if there are previous points project car position to end of time for previous trajectory end
    other_car_s += prev_size*0.02*other_car_speed;

    //if other car is ahead of ego and too close to ego
    // and the front car in target lane is moving slower then ego's speed
    if (other_car_s > car_s && ((other_car_s-car_s)<30 || (other_car_s-car_s)<50 && other_car_speed < speed)) {
      return false;
    }
  }

  if (b_idx != -1) {
    double vx = sensor_fusion[b_idx][3];
    double vy = sensor_fusion[b_idx][4];
    double other_car_speed = sqrt(vx*vx + vy*vy);
    double other_car_s = sensor_fusion[b_idx][5];

    //if other car is behind ego and too close to ego
    // and the back car in target lane is moving faster then ego's speed
    if (other_car_s < car_s && ((car_s-other_car_s)<30 || (car_s-other_car_s)<50 && speed < other_car_speed)) {
      return false;
    }
  }
  return true;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  //variable to hold the current lane. Will start with middle lane i.e. "1"
  int lane = 1;

  //target velocity
  //double ref_vel = 22; //m/s which is equal to 49.5 mph {49.5 / 2.24 = 22}
  double ref_vel = 0.0; //to handle cold start

  h.onMessage([&lane, &ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            // number of data points in previous path still to be consumed by simulator
            int prev_size = previous_path_x.size();

            //the ref is either car_s if there is no previous path
            //or the end of previous_path 
            if (prev_size > 0) {
              car_s = end_path_s;
            }

            bool too_close = false;

            //find if there is a car ahead and it is too close
            for(int i=0; i < sensor_fusion.size(); i++) {

              //car in current lane
              double d = sensor_fusion[i][6];
              if (d >= lane*4. && d <= (lane+1)*4.) {
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double other_car_speed = sqrt(vx*vx + vy*vy);
                double other_car_s = sensor_fusion[i][5];

                //if there are previous points project car position to end of time for previous trajectory end
                other_car_s += prev_size*0.02*other_car_speed;

                //if other car is ahead of ego and too close to ego
                if (other_car_s > car_s && (other_car_s-car_s)<30) {
                  too_close = true;
                  break;
                }
              }
            }

            if (too_close) {
              // if speed is 90% of target, just slow down a bit but do not chnage lane
              // otherwise try to chnage to left lane and otherwise right lane
              if (ref_vel > 0.9*22.0) {
                ref_vel -= 0.1;  //reduce by 0.1 m/s
              } else if (lane > 0 && safetoChangeTo(lane-1, lane, prev_size, car_s, car_speed, sensor_fusion)) {
                //not in leftmost lane and safe to switch to left lane
                lane -= 1;
              } else if (lane < 2 && safetoChangeTo(lane+1, lane, prev_size, car_s, car_speed, sensor_fusion)) {
                //not in rightmost lane and safe to switch to right lane
                lane += 1;
              } else {
                //reduce speed in steps
                ref_vel -= 0.1;  //reduce by 0.1 m/s
              }
            } else if (ref_vel < 22.0) {
              ref_vel += 0.1; //if no car is close in front and ego is below 22 m/s, speed it up
            } 

            // points to crete a spline which will then create actual trajectory points
            vector<double> ptsx;
            vector<double> ptsy;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            //get first two data points either from previous path or from 
            //current car coordinates
            if (prev_size < 2) {

              //get a point which is on tanget of car's path and just before current
              //car position
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              //push two points to vector for spline
              ptsx.push_back(prev_car_x);
              ptsy.push_back(prev_car_y);

              ptsx.push_back(car_x);
              ptsy.push_back(car_y);

            } else {
              //we have atleast two previous point which can be used as first two points of spline
              //let the ref point be the end of previous path
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double ref_x_prev = previous_path_x[prev_size-2];
              double ref_y_prev = previous_path_y[prev_size-2];

              ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

              //push last two points from previous path to spline vector
              ptsx.push_back(ref_x_prev);
              ptsy.push_back(ref_y_prev);
              ptsx.push_back(ref_x);
              ptsy.push_back(ref_y);
            }

            //generate three more points for spline. points are in future
            //generate them in frenet and convert to global X,Y coords
            vector<double> next_wp0 = getXY(car_s+30, 2.0+4.0*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s+60, 2.0+4.0*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s+90, 2.0+4.0*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            //push points to spline vector
            ptsx.push_back(next_wp0[0]);
            ptsy.push_back(next_wp0[1]);
            ptsx.push_back(next_wp1[0]);
            ptsy.push_back(next_wp1[1]);
            ptsx.push_back(next_wp2[0]);
            ptsy.push_back(next_wp2[1]);

            //convert all spline points to ref coordinate system
            //ref coordinate system is either the current car position or
            //last point on prev path

            for (int i = 0; i < ptsx.size(); i++) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = shift_x * cos(ref_yaw) + shift_y * sin(ref_yaw);
              ptsy[i] = - shift_x * sin(ref_yaw) + shift_y * cos(ref_yaw);
            }

            //now create a spline using the 5 points
            tk::spline s;

            //set points in spline
            s.set_points(ptsx, ptsy);

          	//(x,y) points that the car will visit sequentially every .02 seconds
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            //load all previous points into next_x_vals and next_y_vals
            for (int i =0; i < prev_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            //we will now load 50-prev_size points from spline to project the 
            //future path so that we always have 50 points of trajectory being fed

            // the end point of new trajectory in local ref coordinate system
            double target_x = 30; //30 mts straight in local coordinate system
            double target_y = s(target_x);
            double target_dist = distance(0.0, 0.0, target_x, target_y);

            double N = target_dist / (0.02*ref_vel); //each point is visited 0.02 sec 
            double x_delta = target_x/N;
            for(int i =1; i <= 50 - prev_size; i++) {
              //generate (x,y) in local coord
              double x_pt = i * x_delta;
              double y_pt = s(x_pt);

              //convert from local to global coord
              double x_global = x_pt*cos(ref_yaw) -y_pt*sin(ref_yaw) + ref_x;
              double y_global = x_pt*sin(ref_yaw) + y_pt*cos(ref_yaw) + ref_y;

              //push to trajectory vectors
              next_x_vals.push_back(x_global);
              next_y_vals.push_back(y_global);

            }



          	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
