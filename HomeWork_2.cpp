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
    void init_environment(int length, int breadth);
    void init_location_agent(vector<int> agent_location);
    void init_location_goal(vector<int> goal_location);
    double goal_location_length;
    double goal_location_breadth;
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
};



class Simulation{
public:
    void init_simulation(int board_length, int board_breadth,vector<double> agent_location,int case_number,vector<double> goal_location);
    void test_a(vector<double> agent_location, Environment* p_E, Agent* p_A);
    void print_board(Environment* p_E);
    void test_b(Environment* p_E, Agent* p_A);
    void test_c(Environment* p_E, Agent* p_A);
};

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
        if (p_E->goal_location_length > p_A->length_location) {
            //Go down
            p_A->length_location++;
            if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location)=2;
            }
        }else if (p_E->goal_location_breadth < p_A->length_location){
            //Go up
            p_A->length_location--;
            if ((p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) ) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
            }
        }else if (p_E->goal_location_breadth > p_A->breadth_location){
            //Go left
            p_A->breadth_location++;
            if (p_E->board.at(p_A->length_location).at(p_A->breadth_location) != 1) {
                p_E->board.at(p_A->length_location).at(p_A->breadth_location) = 2;
            }
        }else if (p_E->goal_location_breadth < p_A->breadth_location){
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
            
        default:
            break;
    }
    
}




int main(int argc, const char * argv[]) {
    srand(time(NULL));
    int board_length =10;
    int board_breadth =10;
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
    goal_location.push_back(7);
    goal_location.push_back(4);
    Simulation S;
    S.init_simulation(board_length,board_breadth,agent_location,case_number,goal_location);
    return 0;
}
