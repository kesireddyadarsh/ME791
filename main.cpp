//
//  main.cpp
//  ProjectAlpha
//
//  Created by adarsh kesireddy on 2/3/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

#define AKRAND = (double)rand()/RAND_MAX

class Machine{
public:
    double mean;
    double standard_deviation;
    double average;
    void create_machine();
    void pull_arm();
};

//generate single machine
void Machine::create_machine(){
    //generate value for mean
    //generate value for standard deviation
    //Dont do anything with assign values
}

void Machine::pull_arm(){
    //When are is pulled do all the calculation
}

class avl{
public:
    vector<Machine> all_machines;
    vector<Machine>* p_all_machine = &all_machines;
    vector<vector<double>> reward_values;
    vector<vector<double>>* p_reward_values = &reward_values;
    double alpha = 0.1; //Learning rate
    double greedy = 0.25; //greedy level
    
    void create_machines(int number_of_machines);
    void update_averges();
    void leanrning();
};

//generates mutliple machines 
void avl::create_machines(int number_of_machines){
    for (int i=0; i<number_of_machines; i++) {
        Machine M;
        M.create_machine();
        all_machines.push_back(M);
    }
}

void avl::update_averges(){
    //first pull arm
    //second get value from machine
    //Third updare average values of machine
}

void avl::leanrning(){
    int number_of_generations = 100;
    for(int i=0;i<number_of_generations;i++){
        //use update_average function and update averages
        //use learning formule to educate your self
    }
}




int main(int argc, const char * argv[]) {
    srand(time(NULL));
    cout<<"This is start"<<endl;
    return 0;
}
