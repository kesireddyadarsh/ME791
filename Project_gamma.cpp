//
//  main.cpp
//  PorjectGamma
//
//  Created by adarsh kesireddy on 3/24/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <time.h>
#include <stdio.h>

using namespace std;

class City{
public:
    double index;
    double x_position;
    double y_position;
    bool visited;
};

class Environment{
public:
    vector<City> all_cities;
    vector<City>* p_all_cities = &all_cities;
    void create_random_cities(double number_of_cities);
};

void Environment::create_random_cities(double number_of_cities){
    for (int i=0; i<number_of_cities; i++) {
        City C;
        C.index = i;
        //Generate X and Y position
        double temp_x = ((double) rand() / (RAND_MAX))*100;
        double temp_y = ((double) rand() / (RAND_MAX))*100;
        while (temp_x == temp_y) {
            temp_y = rand()%100;
        }
        for (int j=0; j<all_cities.size(); j++) {
            while (temp_x == all_cities.at(j).x_position) {
                temp_x = ((double) rand() / (RAND_MAX))*100;
            }
            while (temp_y == all_cities.at(j).y_position) {
                temp_y = ((double) rand() / (RAND_MAX))*100;
            }
        }
        C.x_position = temp_x;
        C.y_position = temp_y;
        C.visited = false;
        all_cities.push_back(C);
    }
}

class Agent{
public:
    double currentCity;
    vector<double> visited_cities;
    vector<double> distance;
    void calculate_distance_travelled(double next_location, Environment* p_E, Agent* p_A);
};

void Agent::calculate_distance_travelled(double next_location, Environment* p_E, Agent* p_A){
    double dx = p_E->all_cities.at(next_location).x_position - p_E->all_cities.at(p_A->currentCity).x_position;
    double dy = p_E->all_cities.at(next_location).y_position - p_E->all_cities.at(p_A->currentCity).y_position;
    double distance = sqrt((dx*dx)+(dy*dy));
    p_A->distance.push_back(distance);
    cout<<distance<<endl;
}

class Simulation{
public:
    void create_environment();
    Environment E;
    Environment* p_E = &E;
    Agent A;
    Agent* p_A = &A;
    void run_simulation();
};

void Simulation::create_environment(){
    int case_number = 1;
    switch (case_number) {
        case 1:
            p_E->create_random_cities(10);
            break;
            
        case 2:
            p_E->create_random_cities(25);
            break;
            
        case 3:
            p_E->create_random_cities(100);
            break;
    }
    
    cout<<"Print"<<endl;
    for (int i=0; i<p_E->all_cities.size(); i++) {
        cout<<p_E->all_cities.at(i).index<<"\t"<<p_E->all_cities.at(i).x_position<<"\t"<<p_E->all_cities.at(i).y_position<<"\t"<<p_E->all_cities.at(i).visited<<endl;
    }
    run_simulation();
}

void Simulation::run_simulation(){
    //Initialize agent on random city
    //First simulation let it run through all cities. Cities are selected random
    for (int i=0; i<p_E->all_cities.size()-1; i++) {
        //First select random city. change currencity of agent, visited current city to true
        if(i==0){
            p_A->currentCity = rand()%p_E->all_cities.size();
            p_E->all_cities.at(p_A->currentCity).visited = true;
            p_A->visited_cities.push_back(p_A->currentCity);
        }
        //Select next city and make sure its not visited or same as current city
        int next_city = rand()%p_E->all_cities.size();
        while ((next_city == p_A->currentCity) || (p_E->all_cities.at(next_city).visited == true)) {
            next_city = rand()%p_E->all_cities.size();
        }
        //calculate distance travelled
        p_A->calculate_distance_travelled(next_city, p_E, p_A);
        //change current city of agent and visited current city to true
        p_A->currentCity = p_E->all_cities.at(next_city).index;
        p_E->all_cities.at(p_A->currentCity).visited = true;
        p_A->visited_cities.push_back(p_A->currentCity);
    }
    cout<<"Print"<<endl;
    for (int i=0; i<p_E->all_cities.size(); i++) {
        cout<<p_E->all_cities.at(i).index<<"\t"<<p_E->all_cities.at(i).x_position<<"\t"<<p_E->all_cities.at(i).y_position<<"\t"<<p_E->all_cities.at(i).visited<<endl;
    }
    cout<<"Visited cities index"<<endl;
    for (int i=0; i<p_A->visited_cities.size(); i++) {
        cout<<p_A->visited_cities.at(i)<<"\t";
    }
    cout<<endl;
    cout<<"Visited cities index"<<endl;
    for (int i=0; i<p_A->distance.size(); i++) {
        cout<<p_A->distance.at(i)<<"\t";
    }
    cout<<endl;
}

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    cout << "Hello, World!\n";
    Simulation S;
    S.create_environment();
    return 0;
}
