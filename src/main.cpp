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
#include "helper.h"
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

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
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

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
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

vector<double> carCoordinate(double carX, double carY, double carTheta,
							 double worldX, double worldY) {
	vector<double> result(2);
	double d = distance(carX, carY, worldX, worldY);
	double dTheta = atan2(worldY - carY, worldX - carX) - carTheta;

	result[0] = d * cos(dTheta);
	result[1] = d * sin(dTheta);

	return result;
}

vector<vector<double>> splineTrajectoryFromPoints(double car_theta, vector<double> worldX, vector<double> worldY) {
	vector<vector<double>> result(2);

	// transform XY points to car coordinates
	// worldX[1],worldY[1] are car X,Y
	vector<double> X(5);
	vector<double> Y(5);

	for (int i = 0; i < worldX.size(); ++i) {
		// world to car coordinate
		vector<double> rel = carCoordinate(worldX[1], worldY[1], car_theta,
										   worldX[i], worldY[i]);
		X[i] = rel[0];
		Y[i] = rel[1];
	}

	// FIXME PRINTING FUNCTION
	std::cout << "Car points:" << std::endl;
	for (int i = 0; i < X.size(); ++i) {
		std::cout << "(" << X[i] << "," << Y[i] << "), ";
	}
	std::cout << std::endl;

	tk::spline s;
	s.set_points(X, Y);

	// linear approximation
	vector<double> car_points_x;
	vector<double> car_points_y;

	double startX = X[1];
	double startY = Y[1];
	double incX = (X[3] - startX) / 50;  // TODO trajec length variable
	std::cout << "Spline points:" << std::endl;
	for (int i = 0; i < 50; ++i) {  // TODO trajec length variable
		double tmpX = startX + (i + 1) * incX;  // + 1 so first point is after the car, not in the car
		double tmpY = s(tmpX);
		car_points_x.push_back(tmpX);
		car_points_y.push_back(startY + tmpY);
		std::cout << "(" << tmpX << "," << startY + tmpY << ") ";
	}
	std::cout << std::endl;

	// convert car points to world points
	// using carCoordinate function, where car is now world origin in car coordinate system
	vector<double> world_points_x;
	vector<double> world_points_y;
	vector<double> worldOriginFromCar = carCoordinate(worldX[1], worldY[1], car_theta, 0, 0);
	std::cout << "World path points:" << std::endl;
	for (int i = 0; i < car_points_x.size(); ++i) {
		vector<double> tmp = carCoordinate(worldOriginFromCar[0], worldOriginFromCar[1], -car_theta,
										   car_points_x[i], car_points_y[i]);
		world_points_x.push_back(tmp[0]);
		world_points_y.push_back(tmp[1]);
		std::cout << "(" << tmp[0] << "," << tmp[1] << ") ";
	}
	std::cout << std::endl;

	// return path
	result[0] = world_points_x;
	result[1] = world_points_y;
	return result;
}

vector<vector<double>> getTrajectory(int targetLane, double dist_inc, double car_s, double start_theta,
									 double car_x, double car_y,
									 vector<double> map_waypoints_x,
									 vector<double> map_waypoints_y,
									 vector<double> map_waypoints_s) {
	// build point vectors for spline
	// points are:
	//   starting point in past using car heading
	//   starting point
	//   halfway point
	//   end point
	//   end point in future
	vector<double> points_x;
	vector<double> points_y;

	// point in past using heading
	double x0 = car_x - cos(start_theta) * 10 * dist_inc;
	double y0 = car_y - sin(start_theta) * 10 * dist_inc;
	points_x.push_back(x0);
	points_y.push_back(y0);

	// starting point
	points_x.push_back(car_x);
	points_y.push_back(car_y);

	// halfway and endpoint
	double end_s = car_s + 50 * dist_inc;  // TODO make trajectory length variable?
	double end_d = 4 * targetLane + 2;
	vector<double> end_point = getXY(end_s, end_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

	// halfway
	points_x.push_back((car_x + end_point[0]) / 2);
	points_y.push_back((car_y + end_point[1]) / 2);

	// push end point
	points_x.push_back(end_point[0]);
	points_y.push_back(end_point[1]);

	// future end point
	double heading_s = end_s + 10 * dist_inc;  // TODO make heading length variable?
	vector<double> hr = getXY(heading_s, end_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);  // heading reference
	points_x.push_back(hr[0]);
	points_y.push_back(hr[1]);

	// set target values
	double target_x = end_point[0];
	double target_y = end_point[1];

	// theta is zero at east, goes up rotating counter clockwise
	double target_theta = atan2(hr[1] - target_y, hr[0] - target_x);

	// FIXME PRINTING FUNCTION
	std::cout << "Points:" << std::endl;
	for (int i = 0; i < points_x.size(); ++i) {
		std::cout << "(" << points_x[i] << "," << points_y[i] << "), ";
	}
	std::cout << std::endl;

	std::cout << "Starting theta: " << start_theta << ", target theta: " << target_theta << std::endl;

	vector<vector<double>> result = splineTrajectoryFromPoints(start_theta, points_x, points_y);

	return result;
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

  // start in lane 1
  int lane = 1;
  int targetLane = 1;

  // target velocity
  double ref_vel = 49.5;  // mph  TODO use this

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&targetLane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          	double car_yaw_deg = j[1]["yaw"];  // in degrees
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

			// speed limit is 50mph == 22.352 m/s
			double dist_inc = 0.4;  // times 50 is m/s
			double car_yaw = car_yaw_deg * M_PI / 180; // convert to rad

			lane = calculateLane(car_d);

			// vectors to hold the trajectories
			vector<vector<double>> trajectories_next_x_vals;
			vector<vector<double>> trajectories_next_y_vals;

			// vector that hold the target lanes for the generated trajectories
			vector<int> targetLanes;

			// optimal trajectory tracking
			double minCost = 999999999;
			double minIdx;

			// if target lane is not current lane, only append to previous path
			if (lane != targetLane) {
				std::cout << "Moving to a different lane..." << std::endl;
                // TODO increase path size to match 50
                trajectories_next_x_vals.push_back(previous_path_x);
                trajectories_next_y_vals.push_back(previous_path_y);
				minIdx = 0;
				targetLanes.push_back(targetLane);
			} else {
				std::cout << "generating trajectories..." << std::endl;
				// generate up to 3 paths and compare costs: Left, Keep lane, Right
				std::cout << "Keep lane trajectory..." << std::endl;
				vector<vector<double>> klTrajec = getTrajectory(lane, dist_inc, car_s, car_yaw, car_x, car_y,
																map_waypoints_x,
																map_waypoints_y, map_waypoints_s);
				// append to trajectories
				trajectories_next_x_vals.push_back(klTrajec[0]);
				trajectories_next_y_vals.push_back(klTrajec[1]);

				targetLanes.push_back(lane);

				if (lane > 0) { // go left if possible
					std::cout << "Go left trajectory..." << std::endl;
					vector<vector<double>> lTrajec = getTrajectory(lane - 1, dist_inc, car_s, car_yaw, car_x, car_y,
																   map_waypoints_x,
																   map_waypoints_y, map_waypoints_s);
					// append to trajectories
					trajectories_next_x_vals.push_back(lTrajec[0]);
					trajectories_next_y_vals.push_back(lTrajec[1]);

					targetLanes.push_back(lane - 1);
				}

				if (lane < 2) { // go right if possible
					std::cout << "Go right trajectory..." << std::endl;
					vector<vector<double>> rTrajec = getTrajectory(lane + 1, dist_inc, car_s, car_yaw, car_x, car_y,
																   map_waypoints_x,
																   map_waypoints_y, map_waypoints_s);
					// append to trajectories
					trajectories_next_x_vals.push_back(rTrajec[0]);
					trajectories_next_y_vals.push_back(rTrajec[1]);

					targetLanes.push_back(lane + 1);
				}

				// calculate costs
				vector<double> costs;

				for (int i = 0; i < trajectories_next_x_vals.size(); ++i) {
					double cost = calculateCost(car_d, car_s, car_x, car_y, car_yaw, trajectories_next_x_vals[i],
												trajectories_next_y_vals[i], sensor_fusion);
					costs.push_back(cost);
				}

				std::cout << "Costs: ";
				for (int i = 0; i < costs.size(); ++i) {
					std::cout << costs[i] << ", ";

					if (costs[i] < minCost) {
						minCost = costs[i];
						minIdx = i;
					}
				}

				std::cout << std::endl << "Min cost: " << minCost << " at: " << minIdx << std::endl;

			}

			// set target lane to targetLane of picked trajectory
			targetLane = targetLanes[minIdx];

            std::cout << "Sending: lane: " << lane << ", target lane: " << targetLane
					  << ", path length: " << trajectories_next_x_vals[minIdx].size() << "...";
          	msgJson["next_x"] = trajectories_next_x_vals[minIdx];
          	msgJson["next_y"] = trajectories_next_y_vals[minIdx];
            std::cout << "Message sent!" << std::endl;
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
