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

#define population 1000
#define generation 10000
#define cities 100

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
    double startCity;
    vector<double> visited_cities;
    vector<vector<double>> history_visited_cities;
    vector<double> distance;
    vector<double> history_total_distance;
    void calculate_distance_travelled(double next_location, Environment* p_E, Agent* p_A);
    
};


void Agent::calculate_distance_travelled(double next_location, Environment* p_E, Agent* p_A){
    double dx = p_E->all_cities.at(next_location).x_position - p_E->all_cities.at(p_A->currentCity).x_position;
    double dy = p_E->all_cities.at(next_location).y_position - p_E->all_cities.at(p_A->currentCity).y_position;
    double distance = sqrt((dx*dx)+(dy*dy));
    p_A->distance.push_back(distance);
}

class Simulation{
public:
    void create_environment();
    Environment E;
    Environment* p_E = &E;
    Agent A;
    Agent* p_A = &A;
    void init_run_simulation();
    void print_status();
    void evolution();
};

void Simulation::evolution(){
    for (int count_generataions = 0; count_generataions < generation; count_generataions++) {
        //Pick the best routes distance is less
        vector<double> removed_class;
        for (int best_pick_count =0 ; best_pick_count <p_A->history_total_distance.size()/2; best_pick_count++) {
            double rand_1 = rand()%p_A->history_total_distance.size();
            double rand_2 = rand()%p_A->history_total_distance.size();
            while (rand_2 == rand_1) {
                rand_2 = rand()%p_A->history_total_distance.size();
            }
            if (p_A->history_total_distance.at(rand_1) >= p_A->history_total_distance.at(rand_2)) {
                p_A->history_visited_cities.at(rand_1).clear();
                for (int i=0; i<p_A->history_visited_cities.at(rand_2).size(); i++) {
                    p_A->history_visited_cities.at(rand_1).push_back(-100);
                }
                removed_class.push_back(rand_1);
            }else{
                p_A->history_visited_cities.at(rand_2).clear();
                for (int i=0; i<p_A->history_visited_cities.at(rand_1).size(); i++) {
                    p_A->history_visited_cities.at(rand_2).push_back(-100);
                }
                removed_class.push_back(rand_2);
            }
        }
        
        // randomly select a path and place it in empty space. Newly placed path do some random movement
        for (int city_erased = 0; city_erased < removed_class.size(); city_erased++) {
            double rand_1 = rand()%p_A->history_total_distance.size();
            while (p_A->history_visited_cities.at(rand_1).at(0) == -100) {
                rand_1 = rand()%p_A->history_visited_cities.size();
            }
            
            for (int i=0; i<p_A->history_visited_cities.at(rand_1).size(); i++) {
                p_A->history_visited_cities.at(removed_class.at(city_erased)).at(i) = p_A->history_visited_cities.at(rand_1).at(i);
            }
            
            double temp_rand = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
            
            for (int i=0; i<temp_rand; i++) {
                double first_location = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
                double second_location = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
                
                while ((first_location == second_location)||(first_location == 0)||(second_location == 0)) {
                    if (first_location == 0) {
                        first_location = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
                    }
                    if (second_location == 0) {
                        second_location = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
                    }
                    if (first_location == second_location) {
                        first_location = rand()%p_A->history_visited_cities.at(removed_class.at(city_erased)).size();
                    }
                }
                
                double temp_index = p_A->history_visited_cities.at(removed_class.at(city_erased)).at(first_location);
                p_A->history_visited_cities.at(removed_class.at(city_erased)).at(first_location) = p_A->history_visited_cities.at(removed_class.at(city_erased)).at(second_location);
                p_A->history_visited_cities.at(removed_class.at(city_erased)).at(second_location) = temp_index;
            }
            
        }
        
        //Calculate distance
        p_A->history_total_distance.clear();
        removed_class.clear();
        for (int i=0; i<p_A->history_visited_cities.size(); i++) {
            for (int j=0; j<p_A->history_visited_cities.at(i).size()-1; j++) {
                p_A->currentCity = p_A->history_visited_cities.at(i).at(j);
                p_A->calculate_distance_travelled(j+1, p_E, p_A);
            }
            double temp_total_distance=0 ;
            for (int i=0; i<p_A->distance.size(); i++) {
                temp_total_distance += p_A->distance.at(i);
            }
            p_A->history_total_distance.push_back(temp_total_distance);
            p_A->distance.clear();
        }
        
        
        
        cout<<"Distance"<<endl;
        for (int i=0; i<p_A->history_total_distance.size(); i++) {
            cout<<p_A->history_total_distance.at(i)<<"\t";
        }
        cout<<endl;
    }
}

void Simulation::print_status(){
    int debug_mode = 2;
    switch (debug_mode) {
        case 1:
            cout<<"Agent"<<endl;
            for (int i=0; i<p_A->visited_cities.size(); i++) {
                cout<<p_A->visited_cities.at(i)<<"\t";
            }
            cout<<endl;
            for (int i=0; i<p_A->distance.size(); i++) {
                cout<<p_A->distance.at(i)<<"\t";
            }
            cout<<"Environment"<<endl;
            for (int i=0; i<p_E->all_cities.size(); i++) {
                cout<<p_E->all_cities.at(i).index<<"\t"<<p_E->all_cities.at(i).x_position<<"\t"<<p_E->all_cities.at(i).y_position<<"\t"<<p_E->all_cities.at(i).visited<<endl;
            }
            break;
            
        case 2:
            cout<<"Agent"<<endl;
            cout<<"Visited Cities"<<endl;
            for (int i=0; i<p_A->history_visited_cities.size(); i++) {
                for (int j=0; j<p_A->history_visited_cities.at(i).size(); j++) {
                    cout<<p_A->history_visited_cities.at(i).at(j)<<"\t";
                }
                cout<<endl;
            }
            cout<<"Total distance for each route"<<endl;
            for (int i=0; i<p_A->history_total_distance.size(); i++) {
                cout<<p_A->history_total_distance.at(i)<<"\t";
            }
            cout<<endl;
            break;
            
        default:
            break;
    }
    
}

void Simulation::create_environment(){
    p_E->create_random_cities(cities);
    init_run_simulation();
    print_status();
    evolution();
    print_status();
}

void Simulation::init_run_simulation(){
    for (int generate_routes =0; generate_routes < population; generate_routes++) {
        //Initialize agent on random city
        //First simulation let it run through all cities. Cities are selected random
        for (int i=0; i<p_E->all_cities.size()-1; i++) {
            //First select random city. change currencity of agent, visited current city to true
            if(i==0){
                if (generate_routes ==0) {
                    p_A->currentCity = rand()%p_E->all_cities.size();
                    p_A->startCity = p_A->currentCity;
                }
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
        assert(p_A->visited_cities.size() == p_A->distance.size()+1);
        
        p_A->history_visited_cities.push_back(p_A->visited_cities);
        double temp_total_distance=0 ;
        for (int i=0; i<p_A->distance.size(); i++) {
            temp_total_distance += p_A->distance.at(i);
        }
        p_A->history_total_distance.push_back(temp_total_distance);
        p_A->distance.clear();
        p_A->visited_cities.clear();
        p_A->currentCity = p_A->startCity;
        for (int i=0; i<p_E->all_cities.size(); i++) {
            p_E->all_cities.at(i).visited = false;
        }
    }
}

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    cout << "Hello, World!\n";
    Simulation S;
    S.create_environment();
    return 0;
}
