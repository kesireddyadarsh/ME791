//
//  main.cpp
//  ProjectDelta
//
//  Created by adarsh kesireddy on 4/12/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#include "LY_NN.h"

using namespace std;

#define PI 3.14
#define finial_velocity 3.0
#define delta_t 0.2
#define T 5.0
#define x_length_board 1000
#define y_length_board 1000
#define number_of_nn 100


class Goal{
public:
    double goal_1_x_position;
    double goal_1_y_position;
    double goal_2_x_position;
    double goal_2_y_position;
    double goal_mid_point_x;
    double goal_mid_point_y;
    vector<double> box_1;
    vector<double> box_2;
    vector<double> box_3;
    vector<double> box_4;
    vector<double> box_5;
    vector<double> box_6;
    void create_boxes();
};

void Goal::create_boxes(){
    
    //create box_1
    double temp_frent_1_x_position = goal_1_x_position;
    double temp_frent_1_y_position = goal_1_y_position-20;
    double temp_frent_2_x_position = goal_2_x_position;
    double temp_frent_2_y_position = goal_2_y_position-20;
    box_1.push_back(temp_frent_1_x_position);
    box_1.push_back(temp_frent_1_y_position);
    box_1.push_back(goal_1_x_position);
    box_1.push_back(goal_1_y_position);
    box_1.push_back(temp_frent_2_x_position);
    box_1.push_back(temp_frent_2_y_position);
    box_1.push_back(goal_2_x_position);
    box_1.push_back(goal_2_y_position);
    
    //create box_2
    double temp_rear_1_x_position = goal_1_x_position;
    double temp_rear_1_y_position = goal_1_y_position+20;
    double temp_rear_2_x_position = goal_2_x_position;
    double temp_rear_2_y_position = goal_2_y_position+20;
    box_2.push_back(goal_1_x_position);
    box_2.push_back(goal_1_y_position);
    box_2.push_back(temp_rear_1_x_position);
    box_2.push_back(temp_rear_1_y_position);
    box_2.push_back(goal_2_x_position);
    box_2.push_back(goal_2_y_position);
    box_2.push_back(temp_rear_2_x_position);
    box_2.push_back(temp_rear_2_y_position);
    
    //create box_3
    double temp_frent_3_x_position = temp_frent_2_x_position+5;
    double temp_frent_3_y_position = temp_frent_2_y_position;
    double temp_7_x_position = goal_2_x_position+5;
    double temp_7_y_position = goal_2_y_position;
    box_3.push_back(temp_frent_2_x_position);
    box_3.push_back(temp_frent_2_y_position);
    box_3.push_back(goal_2_x_position);
    box_3.push_back(goal_2_y_position);
    box_3.push_back(temp_frent_3_x_position);
    box_3.push_back(temp_frent_3_y_position);
    box_3.push_back(temp_7_x_position);
    box_3.push_back(temp_7_y_position);
    
    //create box_4
    double temp_rear_3_x_position = temp_rear_2_x_position+5;
    double temp_rear_3_y_position = temp_rear_2_y_position;
    box_4.push_back(goal_2_x_position);
    box_4.push_back(goal_2_y_position);
    box_4.push_back(temp_rear_2_x_position);
    box_4.push_back(temp_rear_2_y_position);
    box_4.push_back(temp_7_x_position);
    box_4.push_back(temp_7_y_position);
    box_4.push_back(temp_rear_3_x_position);
    box_4.push_back(temp_rear_3_y_position);
    
    //create box_5
    double temp_frent_4_x_position = temp_rear_1_x_position-5;
    double temp_frent_4_y_position = temp_rear_1_y_position;
    double temp_8_x_position = goal_1_x_position-5;
    double temp_8_y_position = goal_1_y_position;
    box_5.push_back(temp_frent_4_x_position);
    box_5.push_back(temp_frent_4_y_position);
    box_5.push_back(temp_8_x_position);
    box_5.push_back(temp_8_y_position);
    box_5.push_back(temp_frent_1_x_position);
    box_5.push_back(temp_frent_1_y_position);
    box_5.push_back(goal_1_x_position);
    box_5.push_back(goal_1_y_position);
    
    //create box_6
    double temp_rear_4_x_position = temp_frent_1_x_position-5;
    double temp_rear_4_y_positon = temp_rear_1_y_position;
    box_6.push_back(temp_8_x_position);
    box_6.push_back(temp_8_y_position);
    box_6.push_back(temp_rear_4_x_position);
    box_6.push_back(temp_rear_4_y_positon);
    box_6.push_back(goal_1_x_position);
    box_6.push_back(goal_1_y_position);
    box_6.push_back(temp_rear_1_x_position);
    box_6.push_back(temp_rear_1_y_position);
    
    assert(box_1.size() == box_2.size());
    assert(box_1.size() == box_3.size());
    assert(box_1.size() == box_4.size());
    assert(box_1.size() == box_5.size());
    assert(box_1.size() == box_6.size());
    
}

class Ship{
public:
    double ship_x_position;
    double ship_y_position;
    vector<double> x_positions;
    vector<vector<double>> all_paths_x_positions;
    vector<double> y_positions;
    vector<vector<double>> all_paths_y_positions;
    double theta;
    vector<double> theta_values; // Not sure if this will be helpful
    vector<vector<double>> all_thetas_taken;
    double next_speed;
    vector<double> all_speed_from_nn;
    double apple_watch;
    vector<double> fitness;
    bool reached_goal;
    bool box_1;
    bool box_2;
    bool box_3;
    bool box_4;
    bool box_5;
    bool box_6;
};


class Simulation{
public:
    Goal G;
    Goal* p_G = &G;
    Ship S;
    Ship* p_S = &S;
    vector<neural_network> all_networks;
    vector<neural_network>* p_all_networks = &all_networks;
    vector<vector<double>> all_weights;
    vector<vector<double>>* p_all_weights = &all_weights; // for now leave it during EA if not required we can delete it
    void create_environment();
    void run_simulation();
    bool check_if_reached_goal();
    void learning();
    void check_box_number();
    void set_out_of_box();
    void in_box(int box_number);
};

void Simulation::set_out_of_box(){
    S.box_1 = false;
    S.box_2 = false;
    S.box_3 = false;
    S.box_4 = false;
    S.box_5 = false;
    S.box_6 = false;
}

void Simulation::in_box(int box_number){
    switch (box_number) {
        case 1:
            S.box_1 = true;
            break;
        case 2:
            S.box_2 = true;
            break;
        case 3:
            S.box_3 = true;
            S.box_1 = false;
            S.box_2 = false;
            S.box_4 = false;
            S.box_5 = false;
            S.box_6 = false;
            break;
        case 4:
            S.box_4 = true;
            S.box_1 = false;
            S.box_2 = false;
            S.box_3 = false;
            S.box_5 = false;
            S.box_6 = false;
            break;
        case 5:
            S.box_5 = true;
            S.box_1 = false;
            S.box_2 = false;
            S.box_3 = false;
            S.box_4 = false;
            S.box_6 = false;
            break;
        case 6:
            S.box_6 = true;
            S.box_1 = false;
            S.box_2 = false;
            S.box_3 = false;
            S.box_4 = false;
            S.box_5 = false;
            break;
    }
}



void Simulation::check_box_number(){
    if ((S.ship_y_position < G.box_1.at(1))  && (S.ship_x_position < G.box_6.at(0))  ) {
        set_out_of_box();
    }else if ( (S.ship_x_position > G.box_4.at(4))&& (S.ship_y_position > G.box_2.at(3))){
        set_out_of_box();
    }else{
        //Below grid of the goal line
        if (S.ship_y_position >= G.box_1.at(1) && S.ship_y_position <= G.box_1.at(3) ) {
            if ((S.ship_x_position >= G.box_1.at(0) )&&  (S.ship_x_position <= G.box_1.at(4))) {
                in_box(1);
            }else if((S.ship_x_position >= G.box_5.at(0)) && (S.ship_x_position<G.box_5.at(4))){
                in_box(5);
            }else{
                in_box(3);
            }
        }else{ // above grid of goal line
            if ((S.ship_x_position >= G.box_2.at(0)) && (S.ship_x_position <= G.box_2.at(4))) {
                in_box(2);
            }else if ((S.ship_x_position >= G.box_6.at(0))&&(S.ship_x_position<G.box_6.at(4))){
                in_box(6);
            }else{
                in_box(4);
            }
            
        }
    }
}

void Simulation::create_environment(){
    bool testing_mode = true;
    bool dev_mode = false;
    
    //create environment
    
    //create agent location
    if (testing_mode) {
        S.ship_x_position = 200;
        S.ship_y_position = 200;
    }else{
        S.ship_x_position = rand()%x_length_board;
        S.ship_y_position = rand()%y_length_board;
    }
    
    
    //create goal location
    if (testing_mode) {
        G.goal_1_x_position = 500;
        G.goal_1_y_position = 500;
    }else{
        G.goal_1_x_position = rand()%x_length_board;
        while ((G.goal_1_x_position == (S.ship_x_position + 200)) || (G.goal_1_x_position == (S.ship_x_position - 200))) {
            G.goal_1_x_position = rand()%x_length_board;
        }
        G.goal_1_y_position = rand()%y_length_board;
        while ((G.goal_1_y_position == (S.ship_y_position+ 200))||(G.goal_1_y_position == (S.ship_y_position - 200))) {
            G.goal_1_y_position = rand()%y_length_board;
        }
    }
    G.goal_2_x_position = G.goal_1_x_position+5;
    G.goal_2_y_position = G.goal_1_y_position;
    
    if (dev_mode) {
        cout<<"Ship"<<endl;
        cout<<S.ship_x_position<<"\t"<<S.ship_y_position<<endl;
        cout<<"Goal"<<endl;
        cout<<G.goal_1_x_position<<"\t"<<G.goal_1_y_position<<endl;
        cout<<G.goal_2_x_position<<"\t"<<G.goal_2_y_position<<endl;
    }
    
    if (testing_mode) {
        G.goal_mid_point_x = (G.goal_1_x_position+G.goal_2_x_position)/2;
        G.goal_mid_point_y = (G.goal_1_y_position+G.goal_2_y_position)/2;
    }
    
    //create boxes around goal
    G.create_boxes();
    
    //create neuralnetworks with weights
    for (int i=0; i<number_of_nn; i++) {
        neural_network temp_network;
        vector<double> temp_weights;
        temp_network.setup(3, 9, 1);
        for (int j=0; j<46; j++) {
            temp_weights.push_back((double)rand()/RAND_MAX);
        }
        temp_network.set_weights(temp_weights, true);
        p_all_weights->push_back(temp_weights);
        p_all_networks->push_back(temp_network);
    }
    
    assert(p_all_networks->size() == number_of_nn);
    assert(p_all_weights->size() == p_all_networks->size());
    int temp_weight_check = rand()%number_of_nn;
    assert(p_all_weights->at(temp_weight_check).size() == 46);
    
}

bool Simulation::check_if_reached_goal(){
    if (S.box_1 && S.box_2) {
        return true;
    }
    return false;
}

void Simulation::run_simulation(){
    bool dev_mode = true;
    double omega;
    
   //take each network and make it reach goal
    for (int indi_network =0; indi_network<p_all_networks->size(); indi_network++) {
        
        S.reached_goal = false;
        omega = 0;
        
        //set all box check to false
        S.box_1 = false;
        S.box_2 = false;
        S.box_3 = false;
        S.box_4 = false;
        S.box_5 = false;
        S.box_6 = false;
        
        check_box_number();
        
        p_all_networks->at(indi_network).set_in_min_max(0.0, 1000.0);
        p_all_networks->at(indi_network).set_in_min_max(0.0, 1000.0);
        p_all_networks->at(indi_network).set_in_min_max(-180.0, 180.0);
        p_all_networks->at(indi_network).set_out_min_max(-15.0, 15.0);
        
        S.x_positions.push_back(S.ship_x_position);
        S.y_positions.push_back(S.ship_y_position);
        
        //Walk it for number of time steps
        for (int time_step = 0; time_step < 1000; time_step++) {
            //First check if ship is with in board if not kill simulation
            if (S.ship_x_position<0 || S.ship_x_position>x_length_board || S.ship_y_position<0 || S.ship_y_position > y_length_board) {
                //out of board
                S.apple_watch += 1000;
                break;
            }
            
            //Take current x and y position calculate angle
            if (((G.goal_mid_point_y - S.ship_y_position) != 0)||((G.goal_mid_point_x - S.ship_x_position) !=0)) {
                double slope = (G.goal_mid_point_y - S.ship_y_position)/(G.goal_mid_point_x - S.ship_x_position);
                S.theta = atan(slope)*57.29;
            }else if (((G.goal_mid_point_y - S.ship_y_position) == 0)){
                S.theta = 0;
            }else{
                S.theta = 90;
            }
            
            //Pass x,y and angle values into NN and obtain u
            vector<double> input_vector;
            input_vector.push_back(S.ship_x_position);
            input_vector.push_back(S.ship_y_position);
            input_vector.push_back(S.theta);
            p_all_networks->at(indi_network).set_vector_input(input_vector);
            p_all_networks->at(indi_network).execute();
            S.next_speed = p_all_networks->at(indi_network).get_output(0);
            
            //Use equations find new x y positions
            omega += (S.next_speed - omega) *(delta_t/T) ;
            S.theta += omega;
            S.ship_x_position += (finial_velocity * sin(S.theta));
            S.ship_y_position += (finial_velocity * cos(S.theta));
            
            //check if reached goal
            S.reached_goal = check_if_reached_goal();
            S.x_positions.push_back(S.ship_x_position);
            S.y_positions.push_back(S.ship_y_position);
            
            //Check for box
            check_box_number();
            S.reached_goal = check_if_reached_goal();
            
            //Calculate Fitness
            if (!S.reached_goal) {
                S.apple_watch += sqrt((G.goal_mid_point_x - S.ship_x_position)+(G.goal_mid_point_y - S.ship_y_position));
            }else if (S.reached_goal){
                S.apple_watch += -10000;
            }
            
        }
        assert(S.x_positions.size() == S.y_positions.size());
        S.all_paths_x_positions.push_back(S.x_positions);
        S.all_paths_y_positions.push_back(S.y_positions);
        S.fitness.push_back(S.apple_watch);
        cout<<"Network::"<<indi_network<<"\t Fitness::"<<S.apple_watch<<endl;
        S.ship_x_position = S.x_positions.at(0);
        S.ship_y_position = S.y_positions.at(0);
        S.x_positions.clear();
        S.y_positions.clear();
        
    }
    assert(S.all_paths_x_positions.size() == S.all_paths_y_positions.size());
    FILE* print_location;
    print_location=fopen("X and Y locations", "a");
    for (int temp = 0; temp < S.all_paths_x_positions.size(); temp++) {
        for (int temp_1 =0; temp_1 <S.all_paths_x_positions.at(temp).size(); temp_1++) {
            fprintf(print_location, "%f \t %f \n",S.all_paths_x_positions.at(temp).at(temp_1),S.all_paths_y_positions.at(temp).at(temp_1));
        }
        fprintf(print_location, "\n\n");
    }
    fclose(print_location);
    assert(S.fitness.size() == p_all_networks->size());
    
}

void Simulation::learning(){
    //execute EA
    
}


int main(int argc, const char * argv[]) {
    srand(time(NULL));
    Simulation Sim;
    Sim.create_environment();
    Sim.run_simulation();
    return 0;
}
