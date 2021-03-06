using namespace std;

double mph2mstep(double mph_speed) {
    return mph_speed * 0.44704 * 0.02;  // mp/h to m/s to m/step
}

double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int calculateLane(double car_d) {
    /* calculates which lane the car is in given a d coordinate
     *
     * returns -1 when not in a lane
     */

    // 2 is center of left lane, 6 is center of mid lane, 10 is center of right lane
    if (car_d > 1 && car_d < 3) {
        return 0;
    } else if (car_d > 5 && car_d < 7) {
        return 1;
    } else if (car_d > 9 && car_d < 11) {
        return 2;
    } else {
        return -1;
    }
}

double calculateCost(int targetLane, double car_s, double car_d, double car_x, double car_y, double car_theta,
                     vector<double> path_x, vector<double> path_y,
                     vector<vector<double>> other_vehicles) {
    // The data format for each car is: [ id, x, y, vx, vy, s, d]
    double veryCloseW = 6.5;  // per step
    double closeW = 0.8;  // per step
    double slowW = 1.0;
    double laneSwitchW = 1.0;

    // punish driving near other vehicles
    // only consider vehicles in current lane and target lane
    int ownLane = calculateLane(car_d);
    double closeCost = 0.0;
    double veryCloseCost = 0.0;
    float path_mp;  // driving close at the beginning of the trajectory is more dangerous (higher chance of collision)
    for (int i = 0; i < path_x.size(); ++i) {
        path_mp = 1 - i/float(path_x.size());
        for (int vIdx = 0; vIdx < other_vehicles.size(); ++vIdx) {
            int otherVLane = calculateLane(other_vehicles[vIdx][6]);
            if (targetLane == otherVLane || ownLane == otherVLane) {
                double dist = distance(path_x[i], path_y[i], other_vehicles[vIdx][1], other_vehicles[vIdx][2]);
                if (dist < 10) {
                    if (dist < 5.0) {  // one lane width is 4, so should detect cars beside it
                        veryCloseCost += 1.0 * path_mp;
                    } else {
                        closeCost += 1.0 * path_mp;
                    }
                }
            }
        }
    }

    // punish slow driving
    double slowCost = 100 - distance(car_x, car_y, path_x.back(), path_y.back());

    // punish lane switching
    double laneSwitchCost = targetLane == calculateLane(car_d) ? 0.0 : 1.0;

    return veryCloseW * veryCloseCost +
           closeW * closeCost +
           slowW * slowCost +
           laneSwitchW * laneSwitchCost;
}