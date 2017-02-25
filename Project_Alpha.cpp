//
//  main.cpp
//  ProjectAlpha
//
//  Created by adarsh kesireddy on 2/3/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <cstdio>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <assert.h>
#include <math.h>

using namespace std;

const double two_pi = 2.0*3.14159265358979323846;
const double epsilon_box_muller = numeric_limits<double>::min();
const double max_payout_of_any_machine = 100;
const double min_payout_of_any_machine = 0;
const double standard_deviation_of_any_machine = 5;

//create a minimum and maximum payout assign some value
//using those values genereate standard deviation and mean and assign to machine
// reward is equal to pay off
class Machine{
public:
    double mean;
    double standard_deviation;
    double payout_from_machine;
    int number_of_pulls;
    double exepected_payout;
    vector<int> pull_generation;
    
    void create_machine();
    double generate_mean();
    double generate_standard_deviation();
    void reset();
};

double Machine::generate_mean(){
    double temp_rand = (double)rand()/RAND_MAX;
    double temp_mean = (min_payout_of_any_machine+ ( temp_rand * (max_payout_of_any_machine - min_payout_of_any_machine)));
    return temp_mean;
}

double Machine::generate_standard_deviation(){
    double temp_rand = (double)rand()/RAND_MAX;
    double temp_std = temp_rand * standard_deviation_of_any_machine;
    return temp_std;
}

//generate single machine
void Machine::create_machine(){
    //generate value for mean
    mean = generate_mean();
    //generate value for standard deviation
    standard_deviation = generate_standard_deviation();
    reset();
}

void Machine::reset(){
    //Dont do anything with assign values
    number_of_pulls = 0;
    payout_from_machine = 0;
    //temp_reward = 0;
    exepected_payout = max_payout_of_any_machine+20;
}




class MAB{
public:
    vector<Machine> all_machines;
    vector<Machine>* p_all_machine = &all_machines;
    
    void create_machines(int number_of_machines);
};

//generates mutliple machines
void MAB::create_machines(int number_of_machines){
    for (int i=0; i<number_of_machines; i++) {
        Machine M;
        M.create_machine();
        all_machines.push_back(M);
    }
}

class Agent{
public:
    double alpha = 0.1; //Learning rate
    double greedy = 0.25; //greedy level
    int number_of_pulls = 0; // intially number of pulls are zero
    double action_value;
    vector<double> payout_for_graph;
    
    void pull_arm(MAB* p_multiple_bandits);
    void pull_best_arm(MAB* p_multiple_bandits);
    void pull_random_arm(MAB* p_multiple_bandits);
    void calculations(int arm_number,MAB* p_multiple_bandits);
    void learning(int arm_number,MAB* p_multiple_bandits);
    void update_pull_generation(int arm_number,MAB* p_multiple_bandits);
};

void Agent::update_pull_generation(int arm_number,MAB* p_multiple_bandits){
    for (int i=0; i<p_multiple_bandits->p_all_machine->size(); i++) {
        if (i != arm_number) {
            p_multiple_bandits->p_all_machine->at(i).pull_generation.push_back(0);
        }
        if (i == arm_number) {
            p_multiple_bandits->p_all_machine->at(i).pull_generation.push_back(1);
        }
    }
}

void Agent::learning(int arm_number,MAB* p_multiple_bandits){
    p_multiple_bandits->p_all_machine->at(arm_number).exepected_payout = (p_multiple_bandits->p_all_machine->at(arm_number).payout_from_machine*alpha) + (p_multiple_bandits->p_all_machine->at(arm_number).exepected_payout*(1-alpha));
}

void Agent::calculations(int arm_number, MAB* p_multiple_bandits){
    double u1 = 0;
    double u2 = 0;
    double z0 = 0;
    double z1 = 0;
    while (u1 <= epsilon_box_muller)
    {
        u1 = ((double) rand() / (RAND_MAX));
        u2 = ((double) rand() / (RAND_MAX));
    }
    //cout << u1 << endl;
    //cout << u2 << endl;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    //cout << z0 << endl;
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    double sigma = p_multiple_bandits->p_all_machine->at(arm_number).standard_deviation;
    double mu = p_multiple_bandits->p_all_machine->at(arm_number).mean;
    
    p_multiple_bandits->p_all_machine->at(arm_number).payout_from_machine = z0 * sigma + mu;
    payout_for_graph.push_back(p_multiple_bandits->p_all_machine->at(arm_number).payout_from_machine);
}

void Agent::pull_best_arm(MAB* p_multiple_bandits){
    int best_arm;
    double temp_payout;
    for (int i=0; i<p_multiple_bandits->p_all_machine->size(); i++) {
        if (i==0) {
            temp_payout = p_multiple_bandits->p_all_machine->at(i).exepected_payout;
            best_arm = i;
        }else{
            if (temp_payout <= p_multiple_bandits->p_all_machine->at(i).exepected_payout) {
                best_arm = i;
                temp_payout = p_multiple_bandits->p_all_machine->at(i).exepected_payout;
            }
        }
        
    }
    assert(best_arm <= p_multiple_bandits->p_all_machine->size());
    update_pull_generation(best_arm, p_multiple_bandits);
    calculations(best_arm, p_multiple_bandits);
    learning(best_arm, p_multiple_bandits);
}

void Agent::pull_random_arm(MAB* p_multiple_bandits){
    int random_arm = rand()%p_multiple_bandits->p_all_machine->size();
    number_of_pulls++;
    
    p_multiple_bandits->p_all_machine->at(random_arm).number_of_pulls++;
    
    calculations(random_arm, p_multiple_bandits);
    update_pull_generation(random_arm, p_multiple_bandits);
    learning(random_arm, p_multiple_bandits);
}

void Agent::pull_arm(MAB* p_multiple_bandits){
    
    double greedy_random_number = (double)rand()/RAND_MAX;
    
    if (greedy_random_number < (1-greedy)) {
        //work on best arm
        cout<<"Performing on Best arm"<<endl;
        pull_best_arm(p_multiple_bandits);
        
    }else{
        //work on random arm
        cout<<"Perforimg on Random arm"<<endl;
        pull_random_arm(p_multiple_bandits);
    }
    
}

class Simulation{
public:
    void run_simulation(int number_of_arms,int number_of_generations,int case_number);
};

void Simulation::run_simulation(int number_of_arms, int number_of_generations,int case_number){
    
    if (case_number == 1) {
        cout<<"case 1"<<endl;
        MAB mab;
        MAB* p_mab = &mab;
        Agent a;
        Agent* p_a = &a;
        number_of_arms = 1;
        //create simulation
        for (int i=0; i<number_of_arms; i++) {
            Machine M;
            M.create_machine();
            mab.all_machines.push_back(M);
        }
        for (int i=0; i<number_of_generations; i++) {
            if (number_of_arms == 1) {
                p_a->pull_arm(p_mab);
            }else{
                p_a->pull_arm(p_mab);
            }
        }

        assert(p_mab->p_all_machine->at(0).exepected_payout > (p_mab->p_all_machine->at(0).mean - p_mab->p_all_machine->at(0).standard_deviation));
        assert(p_mab->p_all_machine->at(0).exepected_payout < (p_mab->p_all_machine->at(0).mean + p_mab->p_all_machine->at(0).standard_deviation));
        
        cout<<"Test Pass "<<endl;
        
    }else if (case_number == 2)
    {
        cout<<"case 2"<<endl;
        MAB mab;
        MAB* p_mab = &mab;
        Agent a;
        Agent* p_a = &a;
        number_of_arms = 1;
        //create simulation
        for (int i=0; i<number_of_arms; i++) {
            Machine M;
            M.create_machine();
            mab.all_machines.push_back(M);
        }
        for (int i=0; i<number_of_generations; i++) {
            if (number_of_arms == 1) {
                p_a->pull_arm(p_mab);
            }else{
                p_a->pull_arm(p_mab);
            }
        }
        int arm;
        int reward;
        for (int i=0; i<p_mab->p_all_machine->size(); i++ ){
            if (i == 0) {
                reward = p_mab->p_all_machine->at(0).exepected_payout;
                arm = i;
            }else{
                if (reward < p_mab->p_all_machine->at(i).exepected_payout) {
                    reward = p_mab->p_all_machine->at(i).exepected_payout;
                    arm = i;
                }
            }
        }
        for (int i=0; i<p_mab->p_all_machine->size(); i++) {
            assert(p_mab->p_all_machine->at(arm).exepected_payout >= p_mab->p_all_machine->at(i).exepected_payout);
        }
        cout<<"Pass test"<<endl;
    }else{
        cout<<"Case 3:"<<endl;
        MAB multiple_bandits;
        MAB* p_multiple_bandits = &multiple_bandits;
        Agent grand_moff;
        Agent* p_grand_moff = &grand_moff;
        //create simulation
        for (int i=0; i<number_of_arms; i++) {
            Machine M;
            M.create_machine();
            multiple_bandits.all_machines.push_back(M);
        }
        
        for (int stat_run =0; stat_run<30; stat_run++) {
            for (int i=0; i<number_of_generations; i++) {
                if (number_of_arms == 1) {
                    p_grand_moff->pull_arm(p_multiple_bandits);
                }else{
                    p_grand_moff->pull_arm(p_multiple_bandits);
                }
            }
            
            FILE* p_file;
            p_file = fopen("payouts", "a");
            for (int i=0; i<p_grand_moff->payout_for_graph.size(); i++) {
                fprintf(p_file, "%lf \t ",p_grand_moff->payout_for_graph.at(i));
            }
            fprintf(p_file, "\n\n");
            fclose(p_file);
            
            FILE* p_file_2;
            p_file_2 =fopen("GenerationNumber", "a");
            for (int i=0; i<p_multiple_bandits->p_all_machine->size(); i++) {
                for (int j=0; j<p_multiple_bandits->p_all_machine->at(i).pull_generation.size(); j++) {
                    fprintf(p_file_2,"%d \t", p_multiple_bandits->p_all_machine->at(i).pull_generation.at(j) );
                }
                fprintf(p_file_2, "\n\n");
            }
            fprintf(p_file_2, "\n\n\n\n\n");
            fclose(p_file_2);
            
            cout<<"size of payout"<<p_grand_moff->payout_for_graph.size()<<endl;
            for (int rest_number =0; rest_number<p_multiple_bandits->p_all_machine->size(); rest_number++) {
                p_multiple_bandits->p_all_machine->at(rest_number).reset();
                p_multiple_bandits->p_all_machine->at(rest_number).pull_generation.clear();
            }
            p_grand_moff->payout_for_graph.clear();
        }
    }
}

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    cout<<"This is start"<<endl;
    int number_of_arms=3;
    int case_number;
    int number_of_generations = 120;
    cout<<"Enter number of arms"<<endl;
    cin>>number_of_arms;
    cout<<"Enter 1 for Test A\n Enter 2 for Test B\n Enter 3 for Stat run"<<endl;
    cin>>case_number;
    //for learning curve it should be payouofmachine
    
    //cout<<"Number of arms::"<<number_of_arms<<endl;
    Simulation S;
    S.run_simulation(number_of_arms,number_of_generations,case_number);
    return 0;
}
