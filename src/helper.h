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
        std::cout << "Warning, not in lane!";
        return -1;
    }
}

double calculateCost(double car_d, double car_s, double car_x, double car_y, double car_theta,
                     vector<double> path_x, vector<double> path_y,
                     vector<vector<double>> other_vehicles) {
    // The data format for each car is: [ id, x, y, vx, vy, s, d]

    double cost = -1;
    double closePenalty = 0.5;
    double veryClosePenalty = 5.0;

    // punish driving near other vehicles
    for (int i = 0; i < path_x.size(); ++i) {
        for (int vIdx = 0; vIdx < other_vehicles.size(); ++vIdx) {
            if (calculateLane(car_d) == calculateLane(other_vehicles[vIdx][6])) {
                double dist = distance(car_x, car_y, other_vehicles[vIdx][1], other_vehicles[vIdx][2]);
                if (dist < 10) {
                    if (dist < 3) {
                        cost += veryClosePenalty;
                    } else {
                        cost += closePenalty;
                    }
                }
            }
        }
    }

    return cost;
}