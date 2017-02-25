//
//  main.cpp
//  ProjectBeta
//
//  Created by adarsh kesireddy on 2/21/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cassert>
#include <time.h>


using namespace std;

class Environment{
public:
    vector<vector<double>> board;
    vector<int> index;
    void init_environment(int length, int breadth);
    void init_location_agent(vector<int> agent_location);
    void init_location_goal(vector<int> goal_location);
    double goal_location_length;
    double goal_location_breadth;
    vector<vector<double>> Reward;
};

void Environment::init_environment(int length, int breadth){
    vector<double> temp_breadth;
    for (int j=0;j<breadth; j++) {
        temp_breadth.push_back(0);
    }
    for (int i=0; i<length; i++) {
        board.push_back(temp_breadth);
    }
    
    assert(board.size() == length);
    assert(board.at(0).size() == breadth); //Considering its square board
}

void Environment::init_location_agent(vector<int> agent_location){
    board.at(agent_location.at(0)).at(agent_location.at(1)) = 1;
}

void Environment::init_location_goal(vector<int> goal_location){
    board.at(goal_location.at(0)).at(goal_location.at(1)) = 100;
}


class Agent{
public:
    double length_location;
    double breadth_location;
    void move_agent();
    double greedy = 0.1;
    double alpha =0.1;
    double gama = 0.9;
    //Up
    //Down
    //Left
    //Right
    vector<vector<double>> Q_table; //Up Down Left Right
    void sense(Environment* p_E,Agent* p_A);
    void decide(Environment* p_E,Agent* p_A);
    void act(Environment* p_E,Agent* p_A);
    void react(Environment* p_E,Agent* p_A);
    int state_value;
};



class Simulation{
public:
    void init_simulation(int board_length, int board_breadth,vector<double> agent_location,int case_number,vector<double> goal_location);
    void test_a(vector<double> agent_location, Environment* p_E, Agent* p_A);
    void print_board(Environment* p_E);
    void test_b(Environment* p_E, Agent* p_A);
    void test_c(Environment* p_E, Agent* p_A);
    void project(Environment* p_E, Agent* p_A);
    void project_init(Environment* p_E, Agent* p_A);
    void project_run_simulation(Environment* p_E, Agent* p_A);
};

void Simulation::project_run_simulation(Environment* p_E, Agent* p_A){
    double temp_rand = (double)rand()/RAND_MAX;
    if (temp_rand <= (1 - p_A->greedy)) {
        //Do highest value move
        
    }else{
        //Do random value move
        int location = p_A->length_location+(p_A->breadth_location* p_E->board.size());
        //Type of action it needs to perform
        int action = rand()%4;
        //check if action can be performed
        while (p_A->Q_table.at(location).at(action) == -100) {
            action = rand()%4;
        }
        // You know what action to perform on current state
        if (action == 0) {
            //Move up
            p_A->length_location++;
        }else if(action == 1){
            // Move Down
            p_A->length_location--;
        }else if (action == 2){
            //Move Right
            p_A->breadth_location++;
        }else if (action ==3){
            //Move Left
            p_A->breadth_location--;
        }
        //Update location as per action performed
        location = p_A->length_location+(p_A->breadth_location* p_E->board.size());
        //Update Q table
        //Q(s_t, a_t) = Q(s_t, a_t) + alpha_t(r_t+1 + gama*maxQ(s_(t+1), a)- Q(s_t, a_t)
        
    }
}

void Simulation::project_init(Environment* p_E, Agent* p_A){
    vector<double> temp_reward;
    for (int i=0; i < p_E->board.size(); i++) {
        temp_reward.push_back(0);
    }
    for (int i=0; i < p_E->board.size() ; i++) {
        p_E->Reward.push_back(temp_reward);
    }
    
    p_E->Reward.at(p_E->goal_location_length).at(p_E->goal_location_breadth) = 100;
    
    //Inital random numbers for q table
    for (int i=0; i< (p_E->board.size()*p_E->board.size()); i++) {
        vector<double> temp_q_value;
        for (int j=0; j<4; j++) {
            double temp_rand = (((double)rand()/RAND_MAX)*0.1);
            temp_q_value.push_back(temp_rand);
        }
        p_A->Q_table.push_back(temp_q_value);
    }
    
    cout<<"This is Q table"<<endl;
    for (int i=0; i<p_A->Q_table.size(); i++) {
        for (int j=0; j<4; j++) {
            cout<<p_A->Q_table.at(i).at(j)<<"\t";
        }
        cout<<"\n";
    }
    
    //Give values to -100 for movements which agent can't take
    for (int j=0; j<p_E->board.size(); j++) {
        for (int i=0; i<p_E->board.at(j).size(); i++) {
            int location = i+(j*(p_E->board.size()));
            if (j == 0) {
                p_A->Q_table.at(location).at(0) = -100;
            }
            if (j == (p_E->board.size() -1)) {
                p_A->Q_table.at(location).at(1) = -100;
            }
            if (i == (p_E->board.size() -1)) {
                p_A->Q_table.at(location).at(3) = -100;
            }
            if (i==0) {
                p_A->Q_table.at(location).at(2) = -100;
            }
        }
    }
    
    cout<<"This is Q table"<<endl;
    cout<<"Up \t Down \t Left \t Right \n"<<endl;
    for (int i=0; i<p_A->Q_table.size(); i++) {
        for (int j=0; j<4; j++) {
            cout<<p_A->Q_table.at(i).at(j)<<"\t";
        }
        cout<<"\n";
    }
}

void Simulation::project(Environment* p_E, Agent* p_A){
    project_init(p_E, p_A);
    project_run_simulation(p_E,p_A);
}

void Simulation::test_b(Environment* p_E, Agent* p_A){
    while ((p_E->goal_location_breadth != p_A->breadth_location) || (p_E->goal_location_length != p_A->length_location)) {
        int temp_direction;
        cout<<" 1 Right, 2 Left, 3 Up, 4 Down"<<endl;
        cin>>temp_direction;
        switch (temp_direction) {
            case 1:
                p_A->breadth_location++;
                if (p_E->board.at(p_A->length_location).size() <= (p_A->breadth_location)) {
                    p_A->breadth_location--;
                    cout<<"Cant move Right"<<endl;
                }else{
                    if (p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) {
                        p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
                    }
                }
                break;
            case 2:
                p_A->breadth_location--;
                if (0 > p_A->breadth_location) {
                    p_A->breadth_location++;
                    cout<<"Cant move left"<<endl;
                }else{
                    if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                        p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
                    }
                }
                break;
            case 3:
                p_A->length_location--;
                if (0 > p_A->length_location) {
                    p_A->length_location++;
                    cout<<"Cant move Up"<<endl;
                }else{
                    if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                        p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
                    }
                }
                break;
            case 4:
                p_A->length_location++;
                if ( p_E->board.size() <= p_A->length_location) {
                    p_A->length_location--;
                    cout<<"Cant move down"<<endl;
                }else{
                    if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                        p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
                    }
                }
                break;
            default:
                break;
        }
        print_board(p_E);
    }
    cout<<"Position of goal in length and breadth"<<endl;
    cout<<p_E->goal_location_length<<"\t"<<p_E->goal_location_breadth<<endl;
    cout<<"Position of agent in length and breadth"<<endl;
    cout<<p_A->length_location<<"\t"<<p_A->breadth_location<<endl;
    cout<<"Goal Reached!!"<<endl;
}

void Simulation::test_c(Environment* p_E, Agent* p_A){
    while ((p_E->goal_location_breadth != p_A->breadth_location) || (p_E->goal_location_length != p_A->length_location)) {
        if (p_E->goal_location_length >= p_A->length_location) {
            //Go down
            p_A->length_location++;
            if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location)=2;
            }
        }else if (p_E->goal_location_breadth <= p_A->length_location){
            //Go up
            p_A->length_location--;
            if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
            }
        }else if (p_E->goal_location_breadth >= p_A->breadth_location){
            //Go left
            p_A->breadth_location++;
            if (p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
            }
        }else if (p_E->goal_location_breadth <= p_A->breadth_location){
            //Go right
            p_A->breadth_location--;
            if (p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
            }
        }
        print_board(p_E);
    }
    cout<<"Position of goal in length and breadth"<<endl;
    cout<<p_E->goal_location_length<<"\t"<<p_E->goal_location_breadth<<endl;
    cout<<"Position of agent in length and breadth"<<endl;
    cout<<p_A->length_location<<"\t"<<p_A->breadth_location<<endl;
    cout<<"Goal Reached!!"<<endl;
}

void Simulation::print_board(Environment* p_E){
    cout<<"New board look"<<endl;
    for(int i=0;i<p_E->board.size();i++){
        for (int j=0; j<p_E->board.at(i).size(); j++) {
            cout<<p_E->board.at(i).at(j)<<"\t";
        }
        cout<<"\n";
    }
    cout<<endl;
}

void Simulation::test_a(vector<double> agent_location, Environment* p_E, Agent* p_A){
    
    if ((agent_location.at(0) > p_E->board.size()) || (agent_location.at(0) <0)) {
        cout<< "Agent out of board in length. New length position is generated"<<endl;
        int temp = rand()%p_E->board.size();
        agent_location.at(0) = temp;
        cout<<"New location of length is :: \t"<<temp<<endl;
    }
    if ((agent_location.at(1)>p_E->board.at(0).size()) || (agent_location.at(1) <0 )) {
        cout<< "Agent out of board in breadth. New breadth position is generated"<<endl;
        int temp = rand()%p_E->board.at(0).size();
        agent_location.at(1) = temp;
        cout<<"New location of length is :: \t"<<temp<<endl;
    }
    p_E->board.at(agent_location.at(0)).at(agent_location.at(1)) = 1;
    p_A->length_location = agent_location.at(0);
    p_A->breadth_location = agent_location.at(1);
    print_board(p_E);
}

void Simulation::init_simulation(int board_length, int board_breadth,vector<double> agent_location,int case_number,vector<double> goal_location){
    //Initalize board along with placing goal and agent
    Environment E;
    Environment* p_E = &E;
    Agent A;
    Agent* p_A = &A;
    E.init_environment(board_length, board_breadth);
    print_board(p_E);
    p_E->goal_location_length = goal_location.at(0);
    p_E->goal_location_breadth = goal_location.at(1);
    E.board.at(p_E->goal_location_length).at(p_E->goal_location_breadth) = 100;
    switch (case_number) {
        case 1:
            test_a(agent_location, p_E, p_A);
            break;
            
        case 2:
            test_a(agent_location, p_E, p_A);
            test_b(p_E, p_A);
            break;
            
        case 3:
            test_a(agent_location, p_E, p_A);
            test_c(p_E,p_A);
            break;
            
        case 4:
            test_a(agent_location, p_E, p_A);
            project(p_E, p_A);
            break;
            
        default:
            break;
    }
    
}




int main(int argc, const char * argv[]) {
    srand(time(NULL));
    int board_length =5;
    int board_breadth =5;
    //    cout<<"Enter board length"<<endl;
    //    cin>>board_length;
    //    cout<<"Enter board breadth"<<endl;
    //    cin>>board_breadth;
    int case_number ;
    cout<<"Enter 1 to test A\n Enter 2 for Test A and B \n Enter 3 for Test A and C"<<endl;
    cin>>case_number;
    vector<double> agent_location;
    agent_location.push_back(15); //length
    agent_location.push_back(15); //breadth
    vector<double> goal_location;
    goal_location.push_back(2);
    goal_location.push_back(2);
    Simulation S;
    S.init_simulation(board_length,board_breadth,agent_location,case_number,goal_location);
    return 0;
}
