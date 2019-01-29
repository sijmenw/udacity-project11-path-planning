using namespace std;

int calculateLane(double car_d) {
    /* calculates which lane the car is in given a d coordinate
     *
     * returns -1 when not in a lane
     */
    std::cout << std::endl << "car_d: " << car_d << " " << std::endl;
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

double calculateCost(double car_x, double car_y, double car_theta, vector<double> path_x, vector<double> path_y,
                     vector<vector<double>> other_vehicles) {
    // The data format for each car is: [ id, x, y, vx, vy, s, d]

    double cost = -1;
//    cost += rand() % 100;
    // punish strong steering

    // TODO: punish driving near other vehicles

    return cost;
}