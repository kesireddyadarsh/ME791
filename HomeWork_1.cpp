//
//  main.cpp
//  HomeWork_1
//
//  Created by adarsh kesireddy on 1/26/17.
//  Copyright © 2017 AK. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <time.h>
#include <vector>
#include <string>
#include <assert.h>

using namespace std;

// First create a deck which has every card

class card{
public:
    string value;
    string symbol;
    int deck_number;
};

class Deck{
public:
    void create_deck(int decknumber);
    vector<card> allcards;
    vector<card> lefthand;
    vector<card> righthand;
    void shuffledeck();
};

//creates a single deck

void Deck::create_deck(int decknumber){
    //create all the cards
    for (int j= 0 ; j<4; j++) {
        if (j == 0) {
            for (int i =1; i<10; i++) {
                card a;
                a.value = to_string(i+1);
                a.symbol = "hearts";
                allcards.push_back(a);
            }
        }else if (j == 1) {
            for (int i =1; i<10; i++) {
                card a;
                a.value = to_string(i+1);
                a.symbol = "spades";
                allcards.push_back(a);
            }
        }else if (j == 2) {
            for (int i =1; i<10; i++) {
                card a;
                a.value = to_string(i+1);
                a.symbol = "diamonds";
                allcards.push_back(a);
            }
        }else if (j == 3) {
            for (int i =1; i<10; i++) {
                card a;
                a.value = to_string(i+1);
                a.symbol = "clubs";
                allcards.push_back(a);
            }
        }
    }
    
    for (int j=0,temp=0; j<4; j++) {
        for (int i=0; i<4; i++) {
            card a;
            if(i==0){
                a.value="A";
            }else if(i==1){
                a.value="J";
            }else if (i==2){
                a.value="Q";
            }else if (i==3){
                a.value="K";
            }
            if (temp == 0) {
                a.symbol ="hearts";
            }else if (temp == 1){
                a.symbol ="spades";
            }else if (temp == 2){
                a.symbol="diamonds";
            }else if ( temp == 3){
                a.symbol ="clubs";
            }
            allcards.push_back(a);
        }
        temp++;
    }
    
    for (int i=0; i<allcards.size(); i++) {
        allcards.at(i).deck_number = decknumber;
    }
    assert(allcards.size() == 52);
    
}

class shoe{
public:
    void create_shoe(int numberofdecks,int casenumber);
    void shufflecards();
    void bridgeshuffle();
    vector<Deck> allDecks;
    vector<card> shoewithshufflecards;
    vector<card>* pshoewithshufflecards = &shoewithshufflecards;
    vector<Deck>* palldecks = &allDecks;
    vector<card> bridgeshufflecards;
    vector<card> bridgeshufflecards_righthand;
    vector<card> bridgeshufflecards_lefthand;
};

void shoe::bridgeshuffle(){
    for (int i=0; i<shoewithshufflecards.size()/2; i++) {
        bridgeshufflecards_lefthand.push_back(shoewithshufflecards.at(i));
    }
    
    for (int i=(shoewithshufflecards.size()/2); i<shoewithshufflecards.size(); i++) {
        bridgeshufflecards_righthand.push_back(shoewithshufflecards.at(i));
    }
    
    shoewithshufflecards.clear();
    
    assert(bridgeshufflecards_lefthand.size() == bridgeshufflecards_righthand.size());
    
    for (int i=0; i<bridgeshufflecards_righthand.size(); i++) {
        shoewithshufflecards.push_back(bridgeshufflecards_righthand.at(i));
        shoewithshufflecards.push_back(bridgeshufflecards_lefthand.at(i));
    }
    
    bridgeshufflecards_lefthand.clear();
    bridgeshufflecards_righthand.clear();
}

void shoe::create_shoe(int numberofdecks,int casenumber){
    //First generate 4 deck of cards shuffle them
    FILE* pfile;
    pfile = fopen("cards","a");
    for (int i=0; i<numberofdecks; i++) {
        Deck D;
        D.create_deck(i);
        allDecks.push_back(D);
        fprintf(pfile, "Printing Deck number %d without shuffling \n",i);
        for(int j=0; j<allDecks.at(i).allcards.size();j++){
            fprintf(pfile, "%s \t %s\t",palldecks->at(i).allcards.at(j).value.c_str(),palldecks->at(i).allcards.at(j).symbol.c_str());
        }
        fprintf(pfile, "\n");
    }
    
    for (int i=0; i<allDecks.size(); i++) {
        for (int j=0; j<allDecks.at(i).allcards.size(); j++) {
            shoewithshufflecards.push_back(allDecks.at(i).allcards.at(j));
        }
    }
    
    assert(shoewithshufflecards.size() == numberofdecks*52);
    
    switch (casenumber) {
        case 1:
            shufflecards();
            break;
            
        default:
            int loop =0;
            int random_number = rand()%600+1000;
            while (loop<random_number) {
                bridgeshuffle();
                loop++;
            }
            break;
    }
    
    
    for (int i=0; i<shoewithshufflecards.size(); i++) {
        fprintf(pfile, "%s \t %s\t %d\t",shoewithshufflecards.at(i).value.c_str(),shoewithshufflecards.at(i).symbol.c_str(),shoewithshufflecards.at(i).deck_number);
        fprintf(pfile, "\n");
    }
    fclose(pfile);
    
    
}

void shoe::shufflecards(){
    int routate = rand()% 6000 + 10000 ;
    int loop=0;
    //cout<<"THis is routate number::"<<routate<<endl;
    while (loop<routate) {
        //cout<<"This is loop number"<<loop<<endl;
        int rand_1 = rand()%shoewithshufflecards.size();
        int rand_2 = rand()%shoewithshufflecards.size();
        
        while (rand_1 == rand_2) {
            rand_2 = rand()%shoewithshufflecards.size();
        }
        
        swap(shoewithshufflecards.at(rand_2), shoewithshufflecards.at(rand_1));
        loop++;
    }
}

void TESTA(int numberofdecks){
    cout<<" Only one deck is used"<<endl;
    shoe S;
    int casenumber = 2;
    S.create_shoe(numberofdecks,casenumber);
}

void TESTB(){
    int numberofdecks;
    int casenumber;
    cout<<"Enter number of decks:"<<"\t";
    cin>>numberofdecks;
    cout<<"Enter 1 for swap shuffle \n Enter 2 for bridge shuffle \n Enter any other key for bridge shuffle \n"<<"\t";
    cin>>casenumber;
    shoe S;
    S.create_shoe(numberofdecks,casenumber);
}

int main(int argc, const char * argv[]) {
    srand((unsigned)time(NULL));
    
    int testnumber = 2;
    
    if (testnumber == 1) {
        TESTA(testnumber);
    }else if (testnumber ==2 ){
        TESTB();
    }
    
    return 0;
}
