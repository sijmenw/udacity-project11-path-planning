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

vector<vector<double>> splineTrajectoryFromPoints(double car_theta, int startIdx, vector<double> worldX, vector<double> worldY,
												  vector<double> dist_inc_steps) {
	/*
	 * startIdx is the idx of coordinate considered to be the start of the trajectory i.e. 'the car'
	 */
	vector<vector<double>> result(2);
	int trajecL = dist_inc_steps.size();

	// transform XY points to car coordinates
	// worldX[1],worldY[1] are car X,Y
	vector<double> X;
	vector<double> Y;

	for (int i = 0; i < worldX.size(); ++i) {
		// world to car coordinate
		vector<double> rel = carCoordinate(worldX[startIdx], worldY[startIdx], car_theta,
										   worldX[i], worldY[i]);
		X.push_back(rel[0]);
		Y.push_back(rel[1]);
	}

	tk::spline s;
	s.set_points(X, Y);

	// linear approximation using dist_inc_steps as step sizes
	vector<double> car_points_x;
	vector<double> car_points_y;

	double startX = X[startIdx];
	double startY = Y[startIdx];
	double incX = 0;  // the increment, based on step sizes in dist_inc_steps

	for (int i = 0; i < trajecL; ++i) {
		incX += dist_inc_steps[i];
		double tmpX = startX + incX;
		double tmpY = s(tmpX);
		car_points_x.push_back(tmpX);
		car_points_y.push_back(startY + tmpY);
	}

	// convert car points to world points
	// using carCoordinate function, where car is now world origin in car coordinate system
	vector<double> world_points_x;
	vector<double> world_points_y;
	vector<double> worldOriginFromCar = carCoordinate(worldX[startIdx], worldY[startIdx], car_theta, 0, 0);

	for (int i = 0; i < car_points_x.size(); ++i) {
		vector<double> tmp = carCoordinate(worldOriginFromCar[0], worldOriginFromCar[1], -car_theta,
										   car_points_x[i], car_points_y[i]);
		world_points_x.push_back(tmp[0]);
		world_points_y.push_back(tmp[1]);
	}

	// return path
	result[0] = world_points_x;
	result[1] = world_points_y;
	return result;
}

vector<vector<double>> getTrajectory(int targetLane, double dist_inc, double car_s, double start_theta,
									 double car_x, double car_y, double car_speed,
									 vector<double> previous_path_x,
									 vector<double> previous_path_y,
									 vector<double> map_waypoints_x,
									 vector<double> map_waypoints_y,
									 vector<double> map_waypoints_s,
									 int trajecL) {
	// calculate target dist_inc
	// max acc = 10 m/s/s ~ 0.2 m/s/step ~ 0.004 m/step/step
	double max_acc_margin = 0.95;  // to prevent small derivations from exceeding the maximum
	double car_speed_step = mph2mstep(car_speed);
	double max_speed_step = min(dist_inc, car_speed_step + trajecL * 0.004 * max_acc_margin);

	vector<double> dist_inc_steps(trajecL);
	// linearly increase speed from car_speed to max speed
	for (int i = 0; i < trajecL; ++i) {
		// i + 1 so speed increases from first step
		// divided by trajecL because acceleration calculated per step
		dist_inc_steps[i] = car_speed_step + (i+1) * (max_speed_step - car_speed_step) / float(trajecL);
	}

	// build point vectors for spline
	// points are:
	//   starting point in past using car heading
	//   starting point
	// or:
	//   N previous_path points
	//
	// and:
	//   end point
	//   end point in future
	vector<double> points_x;
	vector<double> points_y;
	int buildIdx; // the point considered as the start of the trajectory

	if (previous_path_x.size() < 1) {
		// point in past using heading
		double x0 = car_x - cos(start_theta) * 10 * dist_inc;
		double y0 = car_y - sin(start_theta) * 10 * dist_inc;
		points_x.push_back(x0);
		points_y.push_back(y0);

		// starting point
		points_x.push_back(car_x);
		points_y.push_back(car_y);
		buildIdx = 1;
	} else {
		// assumes previous path is never shorter than 5
		int useNPoints = 5;
		for (int i = 0; i < useNPoints; ++i) {
			points_x.push_back(previous_path_x[i]);
			points_y.push_back(previous_path_y[i]);
		}

		buildIdx = useNPoints - 1;
	}


	// endpoint
	vector<double> end_sd = getFrenet(points_x.back(), points_y.back(), start_theta, map_waypoints_x, map_waypoints_y);
	double end_s = end_sd[0];
	for (int i = 0; i < trajecL; ++i) {
		end_s += dist_inc_steps[i];
	}
	double end_d = 4 * targetLane + 2;
	vector<double> end_point = getXY(end_s, end_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

	// push end point
	points_x.push_back(end_point[0]);
	points_y.push_back(end_point[1]);

	// future end point
	double heading_s = end_s + 10 * dist_inc_steps.back();  // TODO make heading length variable?
	vector<double> hr = getXY(heading_s, end_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);  // heading reference
	points_x.push_back(hr[0]);
	points_y.push_back(hr[1]);

	vector<vector<double>> result = splineTrajectoryFromPoints(start_theta, buildIdx, points_x, points_y, dist_inc_steps);

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
  // speed limit is 50mph
  double ref_vel = 49;  // mph
  double switch_slow = 0.95; // percentage of speed when doing a lane switch
  double speed_reduction_checks = 0.15; // percentage

  int trajecL = 75;
  int trajecMin = 60;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
					  &lane, &targetLane, &ref_vel, &speed_reduction_checks, &switch_slow, &trajecL, &trajecMin]
					  (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

			double dist_inc = mph2mstep(ref_vel);
			double car_yaw = deg2rad(car_yaw_deg);

			int prevLane = lane;
			lane = calculateLane(car_d);
			bool justSwitched;
			if (lane == -1 || prevLane != lane) {
			    justSwitched = true;
			} else {
			    justSwitched = false;
			}

			// vectors to hold the trajectories
			vector<vector<double>> trajectories_next_x_vals;
			vector<vector<double>> trajectories_next_y_vals;
			vector<double> resultX;
			vector<double> resultY;

			// vector that hold the target lanes for the generated trajectories
			vector<int> targetLanes;

			// optimal trajectory tracking
			double minCost = 999999999;
			double minIdx;

            // if path already lengthy enough, only use previous path
			if (previous_path_x.size() >= trajecMin) {
				trajectories_next_x_vals.push_back(previous_path_x);
				trajectories_next_y_vals.push_back(previous_path_y);
				minIdx = 0;
			} else {
				std::cout << "generating trajectories from lane " << lane
                          << " (current target: " << targetLane << ")" << std::endl;
				// generate up to 3 paths and compare costs: Left, Keep lane, Right
				for (int i = 0; i < 3; ++i) {
					vector<vector<double>> klTrajec = getTrajectory(targetLane,
                                                                    (1.0 - i*speed_reduction_checks) * dist_inc,
																	car_s, car_yaw, car_x, car_y, car_speed,
																	previous_path_x, previous_path_y,
																	map_waypoints_x,
																	map_waypoints_y, map_waypoints_s,
																	trajecL);
					// append to trajectories
					trajectories_next_x_vals.push_back(klTrajec[0]);
					trajectories_next_y_vals.push_back(klTrajec[1]);

					targetLanes.push_back(targetLane);
				}

                // only consider lane switching if not already switching
                // allow one iteration for stabilizing
                if (lane == targetLane && justSwitched) {
				    std::cout << "Stabilizing, keeping lane!" << std::endl;
				}
                if (lane == targetLane && !justSwitched) {
                    if (lane > 0) { // go left if possible
                        for (int i = 0; i < 3; ++i) {
                            vector<vector<double>> lTrajec = getTrajectory(lane - 1,
                                                                           switch_slow *
                                                                           (1.0 - i * speed_reduction_checks) *
                                                                           dist_inc,
                                                                           car_s, car_yaw, car_x, car_y, car_speed,
                                                                           previous_path_x, previous_path_y,
                                                                           map_waypoints_x,
                                                                           map_waypoints_y, map_waypoints_s,
                                                                           trajecL);
                            // append to trajectories
                            trajectories_next_x_vals.push_back(lTrajec[0]);
                            trajectories_next_y_vals.push_back(lTrajec[1]);

                            targetLanes.push_back(lane - 1);
                        }
                    }

                    if (lane < 2) { // go right if possible
                        for (int i = 0; i < 3; ++i) {
                            vector<vector<double>> rTrajec = getTrajectory(lane + 1,
                                                                           switch_slow *
                                                                           (1.0 - i * speed_reduction_checks) *
                                                                           dist_inc,
                                                                           car_s, car_yaw, car_x, car_y, car_speed,
                                                                           previous_path_x, previous_path_y,
                                                                           map_waypoints_x,
                                                                           map_waypoints_y, map_waypoints_s,
                                                                           trajecL);
                            // append to trajectories
                            trajectories_next_x_vals.push_back(rTrajec[0]);
                            trajectories_next_y_vals.push_back(rTrajec[1]);

                            targetLanes.push_back(lane + 1);
                        }
                    }
                }

				// calculate costs
				vector<double> costs;

				for (int i = 0; i < trajectories_next_x_vals.size(); ++i) {
					double cost = calculateCost(targetLanes[i], car_s, car_d, car_x, car_y, car_yaw,
												trajectories_next_x_vals[i], trajectories_next_y_vals[i],
												sensor_fusion);
					costs.push_back(cost);
				}

				for (int i = 0; i < costs.size(); ++i) {

					if (costs[i] < minCost) {
						minCost = costs[i];
						minIdx = i;
					}
				}
				std::cout << "Min cost: " << minCost << " at: " << minIdx << std::endl;

				// set target lane to targetLane of picked trajectory
				targetLane = targetLanes[minIdx];

				// start build result trajectory
				if (previous_path_x.size() > 0) {
					for (int i = 0; i < 5; ++i) {
						resultX.push_back(previous_path_x[i]);
						resultY.push_back(previous_path_y[i]);
					}
				}

			}

			// push calculated trajectory onto result trajectory
			for (int i = 0; i < trajectories_next_x_vals[minIdx].size(); ++i) {
				resultX.push_back(trajectories_next_x_vals[minIdx][i]);
				resultY.push_back(trajectories_next_y_vals[minIdx][i]);
			}

          	msgJson["next_x"] = resultX;
          	msgJson["next_y"] = resultY;
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
