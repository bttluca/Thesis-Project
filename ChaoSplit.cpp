
// this program will reproduce the execution of Chao's algorithm modified to take into account the split visit

#define _USE_MATH_DEFINES

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <list>
#include <algorithm>
#include <float.h>

//NOTE:
//1. DOBBIAMO CONSIDERARE IL COSTO DI VISITA SEMPRE SEPARATO DAL COSTO DI VIAGGIO
//2. TODO: DOBBIAMO CONSIDERARE UN LIMITE MASSIMO DI PARTI PER NODO PARI AL NUMERO DI AGENTI

int aaa;
double agentSpeed = 1.0;
double sampleRadius = 1.0;
std::vector<std::vector<double> > travelCosts;
class Node;
class Path;
std::vector<Node*> allNodes;
std::vector<Path*> allPaths;
class LocalDist {
    public:
        double stdDeviation;
        
        //THIS FUNCTION RETURNS THE PROBABILITY TO FIND THE TARGET IN THE GIVEN LOCAL DISTRIBUTION MODEL GIVEN THE AMOUNT OF ENERGY/TIME SPENT IN THE LOCAL SEARCH
        double visitProb(double time) {
            double a = sqrt((agentSpeed * (sampleRadius*2) * time) / M_PI);//DISTANCE FROM CENTRAL POINT OF THE DISTRIBUTION
            return erf(a / (stdDeviation * sqrt(2))); //notice that we don't use the (1/2) * (1+erf(x)) form of the CDF because we consider only the distance from the central point
        }
};

class Partitions {
    public:
        std::vector<std::pair<Path*, std::pair<double, double>>> pathAssignments;//THE FIRST FIELD IS THE PATH, THE SECOND IS THE ENERGY SPENT AND THE EXPECTED SCORE
        Node* node;
};

class Node {
    public:
        int id; //[0-M] value equal to the cell in the cost array that refers to this node
        double score;
        double x; //x coordinate
        double y; //y coordinate
        Partitions parts;//how the node visit is split among paths
        LocalDist model;
        
        Node(int nodeId, double nodeX, double nodeY, double nodeScore) {
            id = nodeId;
            x = nodeX;
            y = nodeY;
            score = nodeScore;
            model = LocalDist();
            model.stdDeviation = 10*sampleRadius;
        }
        
        double visitCost() {
            return (model.stdDeviation*model.stdDeviation*M_PI)/(sampleRadius*sampleRadius*M_PI);
        }
};

class Path {
    public:
        std::list<Node*> nodes;
        double totalScore;
        double energySpent;
        int id;
        
        Path(int pathId) {
            id = pathId;
            totalScore = 0;
            energySpent = 0;
            allPaths.push_back(this);
        }
};

//VARIABILI GLOBALI
double energy;
double deviation;
double record;
int M;
double prob = 5.0;


bool compareTopVsNTop(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.first->id < M);
}

bool comparePartsByEnergySpent(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.second.first > b.second.first);
}

bool comparePartsByEnergyAvailable(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return ((energy - a.first->energySpent + a.second.first) > (energy - b.first->energySpent + b.second.first));
}

bool comparePartsScore(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.second.second > b.second.second);
}

//THIS FUNCTION FIXES THE SCORES OF THE PARTS SO THAT THE SCORES IN TOP PATHS ARE MAXIMIZED AND THE BIGGEST PARTS FAVORED (AS IF THE VISIT OF THOSE WAS EXECUTED FIRST)
void equalizeNode(Node* node) {
    int size = node->parts.pathAssignments.size();
    std::sort(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), compareTopVsNTop);//ORDER PARTS SO THAT TOP PARTS COME FIRST
    int topParts;
    for (int i = 0; i < size; i++) {
        if (node->parts.pathAssignments[i].first->id < M) {
            topParts++;
        }
        else {
            break;
        }
    }
    std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = node->parts.pathAssignments.begin();
    double energyUsed;
    if (topParts > 0) {
        std::advance(it, topParts-1);//THIS ITERATOR POINTS TO THE LAST TOP PATH
        std::sort(node->parts.pathAssignments.begin(), it, comparePartsByEnergySpent);//ORDER TOP PARTS BY ENERGY
    }
    for (int i = 0; i < topParts; i++) {//UPDATE OF THE NTOP PARTS SCORES
        node->parts.pathAssignments[i].second.second = (node->model.visitProb(energyUsed + node->parts.pathAssignments[i].second.first) - node->model.visitProb(energyUsed))*node->score;
        energyUsed = energyUsed + node->parts.pathAssignments[i].second.first;
    }
    it++;//NOW IT POINTS TO THE FIRST NTOP PATH
    std::sort(it, node->parts.pathAssignments.end(), comparePartsByEnergySpent);//ORDER NTOP PARTS BY ENERGY
    for (int i = topParts; i < size; i++) {//UPDATE OF THE NTOP PARTS SCORES
        node->parts.pathAssignments[i].second.second = (node->model.visitProb(energyUsed + node->parts.pathAssignments[i].second.first) - node->model.visitProb(energyUsed))*node->score;
        energyUsed = energyUsed + node->parts.pathAssignments[i].second.first;
    }
}

//TOFIX REFERENCES TO PAIRS
//THIS FUNCTION SPLITS A CERTAIN PARTITION ASSIGNED TO firstPath IN ORDER TO ASSIGN A PART TO secondPath, IT MAXIMIZES THE PART IN secondPath
void splitNode(Node* node, Path* firstPath, Path* secondPath) {
    //TODO: ADD THE CASE WHERE A PART ALREADY EXISTS IN THE TARGET PATH
    std::pair<Path*, std::pair<double, double>> firstPathPart;
    std::pair<Path*, std::pair<double, double>> secondPathPart;
    for (int i = 0; i < node->parts.pathAssignments.size(); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == firstPath) {
            firstPathPart = node->parts.pathAssignments[i];
        }
        if (node->parts.pathAssignments[i].first == secondPath) {
            secondPathPart = node->parts.pathAssignments[i];
        }
    }
    if (secondPathPart.first == NULL) {
        secondPathPart = std::make_pair(secondPath, std::make_pair(0.0, 0.0));
    }
    firstPathPart.second.first = firstPathPart.second.first - std::max(0.0, firstPathPart.second.first - (energy - secondPath->energySpent));
    firstPathPart.first->energySpent = firstPathPart.first->energySpent - std::max(0.0, firstPathPart.second.first - (energy - secondPath->energySpent));
    secondPathPart.second.first = std::min((energy - secondPath->energySpent), firstPathPart.second.first);
    if (firstPathPart.second.first == 0.0) {
        node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), firstPathPart), node->parts.pathAssignments.end());
    }
    equalizeNode(node);       
}

//THIS FUNCTION REVERTS A SPLIT NODE
void revertSplitNode(Node* node, Path* firstPath, Path* secondPath) {
    std::pair<Path*, std::pair<double, double>> firstPathPart;
    std::pair<Path*, std::pair<double, double>> secondPathPart;
    int firstPathFlag = 0;
    int secondPathFlag = 0;
    for (int i = 0; (i < node->parts.pathAssignments.size()) && ((firstPathFlag == 1) && (secondPathFlag == 1)); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == firstPath) {
            firstPathPart = node->parts.pathAssignments[i];
            firstPathFlag = 1;
        }
    }
    for (int i = 0; i < node->parts.pathAssignments.size(); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == secondPath) {
            secondPathPart = node->parts.pathAssignments[i];
            secondPathFlag = 1;
        }
    }
    firstPathPart.second.first = firstPathPart.second.first + secondPathPart.second.first;
    firstPathPart.second.second = firstPathPart.second.second + secondPathPart.second.second;
    node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), secondPathPart), node->parts.pathAssignments.end());
    equalizeNode(node);
}

//this function takes a pointer to a node and returns the cost to reach it and visit it from the start point
double travelCost(Node* a) {
    double x = a->x;
    double y = a->y;
    return sqrt(x*x+y*y)/* + visitCost(a)*/;
}

//this function takes a pointer to a start node "a" and a pointer to a destination node "b" returns the cost to reach it and visit it from the start node
double travelCost(Node* a, Node* b) {
    double x = a->x - b->x;
    double y = a->y - b->y;
    return sqrt(x * x + y * y)/* + visitCost(b)*/;
}

//this function returns the energy cost to visit the node
double visitCost(Node* b) {
    return (b->model.stdDeviation*2)/(sampleRadius*2);
}

//this function takes a vector of pointers to nodes and returns a vector matrix representing the travel+visit costs ordered by node ID number
std::vector<std::vector<double> > updateCosts(std::vector<Node *> allNodes) {
    std::vector<Node *>::iterator it = allNodes.begin();
    //resize the matrix accordingly to the number of nodes
    std::vector<std::vector<double> > distance (allNodes.size());
    for (int i = 0; i < allNodes.size(); i++) {
        distance[i].resize(allNodes.size());
    }
    //update the matrix
    for (int i = 0; i < allNodes.size(); i++) {
        for (int j = 0; j < allNodes.size(); j++) {
            distance[i][j] = travelCost(allNodes[i], allNodes[j]);
            distance[i][j] = distance[i][j];
        } 
    }
    //return statement
    return distance;
}

//this function initializes the nodes by reading their position and scores from a file, it initializes a vector allNodes containing a reference to all the nodes considered.
void initializeNodes(std::string filename) {
    //open file
    std::ifstream myReadFile;
    myReadFile.open(filename);
    //create variables
    std::string line;
    std::string tmp;
    int nodeId = 0;
    double x;
    double y;
    double score;
    Node* newNode;
    if (myReadFile.is_open()) {
        //read line by line
        std::getline(myReadFile, line);
        //divide line in single words
        std::istringstream is( line );
        is >> energy;
        is >> M;
        //first line is the energy available to each agent
        while (!myReadFile.eof()) {
            //the subsequent lines are of the form [x, y, score] where x is the position on the x-axis, y is the position on the y-axis and score is the score associated with the target
            std::getline(myReadFile, line);
            if (myReadFile.eof()) {
                break;
            }
            std::istringstream is( line );
            is >> x;
            is >> y;
            is >> score;
            //create a new node and if its reacheable by the agents add its reference to the allNodes vector, notice that the IDs of nodes are sequential.
            newNode = new Node(nodeId, x, y, score);
            if (travelCost(newNode) < energy) {
                allNodes.push_back(newNode);
                nodeId++;
            }
        }
        //update the travel cost matrix
        travelCosts = updateCosts(allNodes);
    }
    myReadFile.close();
}

//comparison function between nodes to discern which is further from the origin
bool compOrigDistance(Node* a, Node* b) {
    return (sqrt(a->x * a->x + a->y * a->y)) > (sqrt(b->x * b->x + b->y * b->y));
}

//comparison function between paths to discern which has the greater complexive score
bool compPathScore(Path* a, Path* b) {
    return (a->totalScore > b->totalScore);
}


void check() {
    std::cout << "Result:" << std::endl;
    //std::cout << "allPaths.size(): " << allPaths.size() << std::endl;
    for (int i = 0; i < allPaths.size(); i++) {//FOR EACH PATH
        if (allPaths[i]->nodes.size() > 0) {
            std::cout << "path " << allPaths[i]->id << ": ";
            for (std::list<Node*>::iterator it = allPaths[i]->nodes.begin(); it != allPaths[i]->nodes.end(); it++) {//FOR EACH NODE IN THAT PATH
                std::cout << (*it)->id << ": ";
                for (int j = 0; j < (*it)->parts.pathAssignments.size(); j++) {
                    if ((*it)->parts.pathAssignments[j].first == allPaths[i]) {
                        std::cout << "(" << (*it)->parts.pathAssignments[j].second.first << "-" << (*it)->parts.pathAssignments[j].second.second << ")";
                        std::cout << " ";
                        break;
                    }
                    else {
                        std::cout << "rilevato problema: " << (*it)->parts.pathAssignments[j].first->id << "-" << allPaths[i]->id;
                    }
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }
    }
    //std::cin >> aaa;
}

//DONE
//TODO: interrogarsi se non sia possibile avere più posizioni ottime di cui una sia preferibile a causa della possibilità di replanning.
//this function takes a pointer to a node, a pointer to a path and returns a pair composed by the minimal cost of inserting the cost into the path and the optimal insertion position.
std::pair<double, int> fitMove(Node* a, Path* path, double visitCost) {
    double minCost = DBL_MAX;
    double cost = DBL_MAX;
    int position = -1;
    int actual = 0;
    int admissible = 0; //FLAG VARIABLE, 1 == AN ADMISSIBLE SOLUTION WAS ALREADY FOUND, 0 == NO ADMISSIBLE SOLUTION WAS FOUND
    if (path->nodes.size() == 0) {
        cost = travelCost(a)*2 + visitCost;
        if (((cost + path->energySpent) < energy) && ((cost + path->energySpent) < energy)) {//NOTE: cost < minCost IS ALWAYS TRUE IN THIS CASE, SO IT'S OMITTED
            //NOTE THERE'S ONLY ONE POSSIBLE POSITION WHEN path->nodes.size() == 0 SO THE POSITION ASSIGNMENT IS OMITTED
            minCost = cost;
            position = actual;
            admissible = 1;
        }
        return std::make_pair(minCost, 0);
    }
    else {
        //POSITION 0
        std::list<Node*>::iterator it = path->nodes.begin();
        std::list<Node*>::iterator itNext = path->nodes.begin();
        itNext++;
        cost = travelCost(a) + travelCosts[a->id][(*it)->id] - travelCost(*it) + visitCost;
        if (cost + path->energySpent < energy && cost + path->energySpent < energy) {
            minCost = cost;
            position = actual;
            admissible = 1;
        }
        while (itNext != path->nodes.end()) {//POSITION 1 TO LAST
            actual++;
            cost = travelCosts[(*it)->id][a->id] + travelCosts[a->id][(*itNext)->id] - travelCosts[(*it)->id][(*itNext)->id] + visitCost;            
            if (cost < minCost) {
                if (admissible == 0) {//NOTE: IN MOST CASES THE CHECK ON THE ENERGY IS NOT NECESSARY SO WE ADDED THE admissible FLAG TO EXPLOIT LAZY EVALUATION THE NESTED IF IS NECESSARY TO SKIP ASSIGNMENTS IN CASE THE SOLUTION IS UNACCEPTABLE
                    if (cost + path->energySpent < energy) {
                        minCost = cost;
                        position = actual;
                        admissible = 1;
                    }
                }
                else {
                    minCost = cost;
                    position = actual;
                }
            }
            it++;
            itNext++;
        }
        //LAST POSITION
        cost = travelCosts[(*it)->id][a->id] + travelCost(a) - travelCost(*it);
        if (cost < minCost) {
            if (admissible == 0) {//NOTE: IN MOST CASES THE CHECK ON THE ENERGY IS NOT NECESSARY SO WE ADDED THE admissible FLAG CHECK TO SPARE SOEM CYCLES
                if (cost + path->energySpent < energy) {
                    minCost = cost;
                    position = actual;
                    //admissible = 1; //NOTE: IT'S NONSENSE WASTING TIME TO SET THE FLAG WHEN WE ALREADY REACHED THE END OF THE PATH
                }
            }
            else {
                minCost = cost;
                position = actual;
            }
        }
        return std::make_pair(minCost, position);
    }
}

//this function takes a vector of nodes, a vector of paths and returns the cheapest insertion of a node from the first vector to a path in the second in the form of pair of pairs. The first pair is the chosen node's id and the cost. The second one is the path and the position where to insert.
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = -1;
    int next = -1;
    double cost;
    int admissible;
    int samePath;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//CHECK ON EACH POSITION
            admissible = 0;
            samePath = 0;
            for (int l = 0; l < nodes[i]->parts.pathAssignments.size(); l++) {//FOR EACH PART BELONGING TO THE NODE 
                if ((nodes[i]->parts.pathAssignments[l].first = paths[j]) && (nodes[i]->parts.pathAssignments[l].first != NULL)) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {
                std::pair<double, int> result = fitMove(nodes[i], paths[j], visitCost(nodes[i]));//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if (cost < minCost && paths[j]->energySpent + cost < energy) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths[j]).id;
                        minNode = nodes[i]->id;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths[j]->id;
                        minNode = nodes[i]->id;
                    }
                }
            }
        }
    }
    return std::make_pair(std::make_pair(minNode, minCost), std::make_pair(minPath, position));
}

//this function takes a vector of nodes, a vector of paths and returns the cheapest insertion of a node from the first vector to a path in the second in the form of pair of pairs. The first pair is the chosen node's id and the cost. The second one is the path and the position where to insert.
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths, std::pair<Path*, std::pair<double, double>>* part) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = -1;
    int next = -1;
    int samePath = 0;;
    double cost;
    int admissible;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//CHECK ON EACH POSITION
            admissible = 0;
            samePath = 0;
            for (std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = nodes[i]->parts.pathAssignments.begin(); it != nodes[i]->parts.pathAssignments.end(); it++) {
                if ((*it).first = paths[j]) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {
                std::pair<double, int> result = fitMove(nodes[i], paths[j], part->second.first);//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if (cost < minCost && paths[j]->energySpent + cost < energy) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths[j]).id;
                        minNode = i;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths[j]->id;
                        minNode = nodes[i]->id;
                    }
                }
            }
        }
    }
    return std::make_pair(std::make_pair(minNode, minCost), std::make_pair(minPath, position));
}

//THIS FUNCTION TAKES A REFERENCE TO THE PART MOVED, AN INTEGER REPRESENTING THE ID OF THE NODE, AN INT REPRESENTING THE ID OF THE DESTINATION PATH AND AN INT WHICH IS THE POSITION IN THE NEW PATH. THEN IT MOVES THE NODE WITH SUCH AN ID IN THE POSITION IN THE NEW PATH AND UPDATES THE path FIELD IN THE NODE
void movePart(std::pair<Path*, std::pair<double, double>>* part, int nodeId, int destPath, int position) {
    if (part->first != NULL) {//IF THE PREVIOUS PATH IS NOT NULL REMOVE PART FROM THE PREVIOUS PATH
        Path* originPath = part->first;
        part->first = NULL;
        std::list<Node*>::iterator it = find(originPath->nodes.begin(), originPath->nodes.end(), allNodes[nodeId]);
        originPath->totalScore = originPath->totalScore - part->second.second;
            std::list<Node*>::iterator itNext = it;
            itNext++;
        if (it != originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCosts[(*itPrev)->id][nodeId] - travelCosts[nodeId][(*itNext)->id] + travelCosts[(*itPrev)->id][(*itNext)->id];
        }
        if (it == originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCost(allNodes[nodeId]) - travelCosts[nodeId][(*itNext)->id];
        }
        if (it != originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            originPath->energySpent = originPath->energySpent - travelCosts[(*itPrev)->id][nodeId] - travelCost(allNodes[nodeId]);
        }
        if (it == originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            originPath->energySpent = originPath->energySpent - (travelCost(allNodes[nodeId])*2);
        }
        originPath->nodes.erase(std::remove(originPath->nodes.begin(), originPath->nodes.end(), allNodes[nodeId]), originPath->nodes.end());
    }
    //ADD PART TO THE NEW PATH
    if (destPath >= 0) {//IF THE id OF THE DESTINATION PATH IS POSITIVE THE DESTINATION PATH IS ASSUMED TO EXIST
        std::list<Node*>::iterator it = allPaths[destPath]->nodes.begin();
        std::advance(it, position);
        part->first = allPaths[destPath];
        allPaths[destPath]->nodes.insert(it, allNodes[nodeId]);
        it = allPaths[destPath]->nodes.begin();
        std::advance(it, position);
        allPaths[destPath]->totalScore = allPaths[destPath]->totalScore + part->second.second;
        if (position != 0 && position != allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCosts[(*itPrev)->id][nodeId] + travelCosts[nodeId][(*itNext)->id];
        }
        if (position == 0 && position != allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCost(allNodes[nodeId]) + travelCosts[nodeId][(*itNext)->id];
        }
        if (position != 0 && position == allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCosts[(*itPrev)->id][nodeId] + travelCost(allNodes[nodeId]);
        }
        if (position == 0 && position == allPaths[destPath]->nodes.size()-1) {
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCost(allNodes[nodeId])*2;
        }
    }
}

//THIS FUNCTION CONSIDERS THE PARTITIONS OF A NODE AND TRIES TO MERGE THE PARTS IN A NTOP PATH TO THE PARTS IN A TOP PATH
void mergeNode(Node* node) {
    int size = node->parts.pathAssignments.size();
    std::sort(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), compareTopVsNTop);//ORDER PARTS SO THAT TOP PARTS COME FIRST
    int topParts = 0;
    for (int i = 0; i < size; i++) {
        if (node->parts.pathAssignments[i].first->id < M) {
            topParts++;
        }
        else {
            break;
        }
    }
    std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = node->parts.pathAssignments.begin();
    if (topParts > 0) {
        std::advance(it, topParts-1);// std::(it, topParts-1);//THIS ITERATOR POINTS TO THE LAST TOP PART
        std::sort(node->parts.pathAssignments.begin(), it, comparePartsByEnergyAvailable);//NOTE: ORDER TOP PARTS BY SUM OF ENERGY LEFT IN THE PATH AND ENERGY ALREADY ASSIGNED TO THE VISIT
        it++;//NOW IT POINTS TO THE FIRST NTOP PATH
    }
    std::sort(it, node->parts.pathAssignments.end(), comparePartsByEnergyAvailable);//NOTE: ORDER NTOP PARTS BY SUM OF ENERGY LEFT IN THE PATH AND ENERGY ALREADY ASSIGNED TO THE VISIT
    for (int i = 0; i < topParts; i++) {
        for (int j = size - 1; j > (topParts - 1); j--) {
            node->parts.pathAssignments[i].second.first = node->parts.pathAssignments[i].second.first + std::min((energy - node->parts.pathAssignments[i].first->energySpent), node->parts.pathAssignments[j].second.first);
            node->parts.pathAssignments[j].second.first = std::max(0.0, (node->parts.pathAssignments[j].second.first - (energy - node->parts.pathAssignments[i].first->energySpent)));
            if ((energy - node->parts.pathAssignments[i].first->energySpent) == 0) {
                break;
            }
        }
        for (int j = topParts - 1; j > i; j--) {
            node->parts.pathAssignments[i].second.first = node->parts.pathAssignments[i].second.first + std::min((energy - node->parts.pathAssignments[i].first->energySpent), node->parts.pathAssignments[j].second.first);
            node->parts.pathAssignments[j].second.first = std::max(0.0, (node->parts.pathAssignments[j].second.first - (energy - node->parts.pathAssignments[i].first->energySpent)));
            if ((energy - node->parts.pathAssignments[i].first->energySpent) == 0) {
                break;
            }
        }
    }
    for (int i = 0; i < topParts; i++) {
        if (node->parts.pathAssignments[i].second.first == 0.0) {
            movePart(&node->parts.pathAssignments[i], node->id, -1, 0);
            node->parts.pathAssignments[i].first->nodes.erase(std::remove(node->parts.pathAssignments[i].first->nodes.begin(), node->parts.pathAssignments[i].first->nodes.end(), node), node->parts.pathAssignments[i].first->nodes.end());
            node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), node->parts.pathAssignments[i]), node->parts.pathAssignments.end());
        }
    }
    equalizeNode(node);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, std::vector<Path*> paths, std::pair<Path*, std::pair<double, double>>* part) {
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths, part);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, std::vector<Path*> paths) {
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, Path* path) {
    std::vector<Path*> paths;
    paths.push_back(path);
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(std::vector<Node*> nodes, Path* path) {
    std::vector<Path*> paths;
    paths.push_back(path);
    return cheapestInsertion(nodes, paths);
}

//BOOLEAN FUNCTION TO SORT PATHS BY SCORE (GREATER SCORE FIRST)
bool scorePathComp(Path* a, Path* b) {
    return a->totalScore > b->totalScore;
}

//THIS FUNCTION GETS A COPY OF THE topPaths VECTOR AND SUMS THE SCORES OF THE PATHS IN IT
double calculateRecord(std::vector<Path*> paths){
    std::vector<Path*> tmpAllPaths = paths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    double sum = 0.0;
    for (int i = 0; i < M; i++) {
        sum = sum + tmpAllPaths[i]->totalScore;
    }
    return sum;
}

//THIS FUNCTION CREATES THE INITIAL PATHS. IT ADDS THE FURTHEST NODES TO A FEW PATHS AND THEN ADDS THE REMAINING NODES WITH CHEAPEST INSERTION
void initializePaths(std::vector<Node*> nodesToAdd) {
    int nodesNumber = nodesToAdd.size();
    int L = std::min(5, nodesNumber);//NOTE: IN THE ORIGINAL WORK PRESENTED BY CHAO THE NUMBER OF AGENTS IS LOWER THAN 5. IF WE GO ABOVE WE MAY NEED TO CHANGE THIS LINE OF CODE.
    Path* pathPtr;
    //sort the nodes in descending order of distance from origin
    std::sort(nodesToAdd.begin(), nodesToAdd.end(), compOrigDistance);
    for (int i = 0; i < M && i < L; i++) {//FOR EACH INITIAL PATH (THE NUMBER OF INITIAL PATH IS THE NUMBER OF AGENTS) ADD THE FURTHEST NODE FROM THE START
        pathPtr = new Path(i);
        nodesToAdd[i]->parts.pathAssignments.push_back(std::make_pair(pathPtr, std::make_pair(visitCost(nodesToAdd[i]), nodesToAdd[i]->score)));
        movePart(&nodesToAdd[i]->parts.pathAssignments[0], nodesToAdd[i]->id, pathPtr->id, pathPtr->nodes.size());
        nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), nodesToAdd[i]), nodesToAdd.end());
    }
    int prevSize;
//INFINITE LOOP BELOW, APPARENTLY THE nodesToAdd SIZE DIMINISHES ONLY IN THE FIRST CYCLE
    while (nodesToAdd.size() > 0) {//WHILE THERE ARE NODES WITHOUT PATH WE INSERT ONE OF THOSE IN A PATH, WHEN NECESSARY WE CREATE A NEW ONE
        prevSize = nodesToAdd.size();
        std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(nodesToAdd, allPaths);
        if (result.first.first == -1) {//IF THERE'S NO ADMISSIBLE INSERTION THEN CREATE A NEW PATH AND REPEAT
            pathPtr = new Path(allPaths.size());
        }
        else {//ELSE ADD THE NODE TO THE PATH, UPDATE THE NODE'S PATH FIELD AND THE PATH'S SCORE AND ENERGY SPENT FIELDS
            allNodes[result.first.first]->parts.pathAssignments.push_back(std::make_pair(allPaths[result.second.first], std::make_pair(visitCost(allNodes[result.first.first]), allNodes[result.first.first]->score)));
            movePart(&allNodes[result.first.first]->parts.pathAssignments[0], result.first.first, result.second.first, result.second.second);
            nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), allNodes[result.first.first]), nodesToAdd.end());
        }
    }
    record = calculateRecord(allPaths);//SOME NON-TRIVIAL WORK WAS DONE ON THE allPaths VECTOR SO WE UPDATE THE record
}

//ARRIVATO QUI
//NOTE:
//1. DOBBIAMO CONSIDERARE IL COSTO DI VISITA SEMPRE SEPARATO DAL COSTO DI VIAGGIO
//2. QUANDO DECIDIAMO SE SPLITTARE UN NODO DOBBIAMO CONSIDERARE SE NE VALE LA PENA, SUGGERISCO DI IMPOSTARE UNA PERCENTUALE PARI A QUELLA DELLA DEVIATION E SE UNA QUALCHE FORMA DI RAPPORTO SCORE/COSTO SCENDE SOTTO QUELLA PERCENTUALE NON ESEGUIRE LO SPLIT
//TODO (MAYBE): ADD A METHOD THAT UPDATES THE SCORE OF PATHS
//TODO: ADD A GLOBAL VECTOR TOP PATHS AND A GLOBAL VECTOR NTOP PATHS SO THAT THEY'RE ALWAYS ACCOUNTED FOR CORRECTLY
//TODO: ADD A isTopPath(Path* path) BOOLEAN FUNCTION THAT CHECKS IF path IS A TOP PATH AND USE IT INSTEAD OF CHECKING THE PATH ID TO KNOW IF THE PATH IS A TOP PATH
//TODO: FIX THE WHOLE allPaths[] ORDER SITUATION BY LETTING THE calculateRecord() CREATE A COPY OF THE VECTOR TO ORDER AND CALCULATE THE RECORD ON THAT. OR BY ADDING A topPaths and a nTopPaths VECTOR.
//TODO: AGGIUNGERE UN CASO A splitNode NEL CASO IN CUI IL NODO DI CUI FACCIAMO LO SPLIT HA GIA' UNA PARTE NEL DESTINATION PATH

//THIS FUNCTION IMPLEMENTS THE ONE-POINT MOVEMENT AS OUTLINED IN THE LITERATURE
void onePointMove() {
    int bestMovePosition;
    int bestMovePath = -1;
    std::pair<Path*, std::pair<double, double>>* bestMovePart;
    std::pair<Path*, std::pair<double, double>>* candPart;
    double bestMove = 0;
    double cost = 0;
    int actionDone;
    int prevPath;
    int prevPosition;
    double candRecord;
    double totScore; 
    for (int i = 0; i < allNodes.size(); i++) {//FOR EACH NODE
        bestMovePath = -1;
        bestMovePosition = -1;
        candRecord = -DBL_MAX;
        actionDone = 0;
        mergeNode(allNodes[i]);
        for (int l = 0; (l < allNodes[i]->parts.pathAssignments.size()) && (actionDone == 0); l++) {//FOR EACH PART
            candPart = &allNodes[i]->parts.pathAssignments[l];
            prevPath = candPart->first->id;
            prevPosition = 0;
            for (std::list<Node*>::iterator tmpIt = allPaths[prevPath]->nodes.begin(); tmpIt != allPaths[prevPath]->nodes.end(); tmpIt++) {
                if (*(tmpIt) == allNodes[i]) {//IF WE'VE FOUND THE POSITION WE WXIT THE FOR LOOP
                    break;
                }
                prevPosition++;
            }
            for (int j = 0; (j < allPaths.size()) && (actionDone == 0); j++) {//WE TRY TO MOVE IT IN ANOTHER PATH
                Path* candPath = allPaths[j];
                if (candPath->id != prevPath) {
                    std::pair<double, int> result = fitMove(allNodes[i], candPath, candPart->second.first);
                    cost = result.first;
                    if (cost + candPath->energySpent < energy) {
                        movePart(candPart, i, j, result.second);
                        totScore = calculateRecord(allPaths);
                        if (totScore > record) {
                            record = totScore;
                            deviation = (1.0-(prob/100))*record;
                            actionDone = 1;
                        }
                        else {
                        movePart(candPart, i, prevPath, prevPosition);
                            if (totScore > deviation && totScore > candRecord) {//IF THE SOLUTION IS NOT A STRICT IMPROVEMENT WE SAVE THE BEST SOLUTION THAT DOES NOT WORSEN THE SOLUTION BY MORE THAN THE DEVIATION
                                bestMovePath = j;
                                bestMovePosition = result.second;
                                candRecord = totScore;
                                bestMovePart = &allNodes[i]->parts.pathAssignments[l];
                            }
                        }
                    }
                    else {
                        //SAVE PATH AND NODE STATUS
                        Node tmpNode = *(allNodes[i]);
                        Path tmpPath1 = *(allPaths[j]);
                        Path tmpPath2 = *(allPaths[prevPath]);
                        //PERFORM SPLIT NODE
                        if ((result.first / candPart->second.second) > (prob/100)) {
                            splitNode(allNodes[i], allPaths[prevPath], allPaths[j]);
                            totScore = calculateRecord(allPaths);
                        }
                        if (totScore > record) {//CONFIRM SPLIT NODE
                            record = totScore;
                            deviation = (1.0-(prob/100))*record;
                            actionDone = 1;
                        }
                        else {//REVERT SPLIT NODE
                            *(allNodes[i]) = tmpNode;
                            *(allPaths[j]) = tmpPath1;
                            *(allPaths[prevPath]) = tmpPath2;
                        }
                    }
                }
            }
        }
        if (bestMovePath != -1 && actionDone == 0) {
            movePart(bestMovePart, i, bestMovePath, bestMovePosition);
            record = calculateRecord(allPaths);
            deviation = (1.0-(prob/100))*record;
        }
    }
}

//COMPARISON FUNCTION BETWEEN NODES BASED ON SCORE
bool compareNodeScore(Node* a, Node* b) {
    if (a->score < b->score) {
        return false;
    }
    else {
        return true;
    }
}

//COMPARISON FUNCTION BETWEEN PARTS BASED ON SCORE
/*
bool comparePartScore(std::pair<Path*, std::pair<double, double>>* a, std::pair<Path*, std::pair<double, double>>* b) {
    if (a->second.second < b->second.second) {
        return false;
    }
    else {
        return true;
    }
}*/

bool comparePartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>*> a, std::pair<int, std::pair<Path*, std::pair<double, double>>*> b) {
    if (a.second->second.second < b.second->second.second) {
        return false;
    }
    else {
        return true;
    }
}

//COMPARISON FUNCTION BETWEEN NODES BASED ON SCORE WEIGHTED ON THE LOCAL COST OF THAT NODE
/*bool compareWeightedNodeScore(Node* a, Node* b) {
    std::list<Node*>::iterator itA =  find(a->path->nodes.begin(), a->path->nodes.end(), a);
    std::list<Node*>::iterator aNext = itA;
    std::list<Node*>::iterator aPrev = itA;
    std::list<Node*>::iterator itB =  find(b->path->nodes.begin(), b->path->nodes.end(), b);
    std::list<Node*>::iterator bNext = itA;
    std::list<Node*>::iterator bPrev = itA;
    aNext++;
    aPrev--;
    bNext++;
    bPrev--;
    double costA = travelCosts[(*aPrev)->id][a->id] + travelCosts[a->id][(*aNext)->id] - travelCosts[(*aPrev)->id][(*aNext)->id];
    double costB = travelCosts[(*bPrev)->id][a->id] + travelCosts[a->id][(*bNext)->id] - travelCosts[(*bPrev)->id][(*bNext)->id];
    return (a->score/costA < b->score/costB);
}*/


//COMPARISON FUNCTION BETWEEN PARTS BASED ON SCORE WEIGHTED ON THE LOCAL COST OF THAT PART
bool compareWeightedPartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>*> a, std::pair<int, std::pair<Path*, std::pair<double, double>>*> b) {
    std::list<Node*>::iterator itA =  find(a.second->first->nodes.begin(), a.second->first->nodes.end(), allNodes[a.first]);
    std::list<Node*>::iterator itALast = a.second->first->nodes.end();
    itALast--;
    std::list<Node*>::iterator itB =  find(b.second->first->nodes.begin(), b.second->first->nodes.end(), allNodes[b.first]);
    std::list<Node*>::iterator itBLast = b.second->first->nodes.end();
    itBLast--;
    double costA = 0.0;
    if ((itA != a.second->first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        aNext++;
        double costA = travelCosts[(*aPrev)->id][(*itA)->id] + travelCosts[(*itA)->id][(*aNext)->id] - travelCosts[(*aPrev)->id][(*aNext)->id];
    }
    if ((itA != a.second->first->nodes.begin()) && (itA == itALast)) {
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        double costA = travelCosts[(*aPrev)->id][(*itA)->id] + travelCost(allNodes[(*itA)->id]) - travelCost(allNodes[(*aPrev)->id]);
    }
    if ((itA == a.second->first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        aNext++;
        double costA = travelCost(allNodes[(*itA)->id]) + travelCosts[(*itA)->id][(*aNext)->id] - travelCost(allNodes[(*aNext)->id]);
    }
    if ((itA == a.second->first->nodes.begin()) && (itA == itALast)) {
        double costA = travelCost(allNodes[(*itA)->id]) + travelCost(allNodes[(*itA)->id]);
    }
    double costB = 0.0;
    if ((itB != b.second->first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        bNext++;
        double costB = travelCosts[(*bPrev)->id][(*itB)->id] + travelCosts[(*itB)->id][(*bNext)->id] - travelCosts[(*bPrev)->id][(*bNext)->id];
    }
    if ((itB != b.second->first->nodes.begin()) && (itB == itBLast)) {
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        double costB = travelCosts[(*bPrev)->id][(*itB)->id] + travelCost(allNodes[(*itB)->id]) - travelCost(allNodes[(*bPrev)->id]);
    }
    if ((itB == b.second->first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        bNext++;
        double costB = travelCost(allNodes[(*itB)->id]) + travelCosts[(*itB)->id][(*bNext)->id] - travelCost(allNodes[(*bNext)->id]);
    }
    if ((itB == b.second->first->nodes.begin()) && (itB == itBLast)) {
        double costB = travelCost(allNodes[(*itB)->id]) + travelCost(allNodes[(*itB)->id]);
    }
    return (a.second->second.second/costA < b.second->second.second/costB);
}

//THIS FUNCTION IMPLEMENTS THE CLEAN UP ROUTINE
void cleanUp() {
    double record = 0.0;
    int actionDone = 0;
    double cost = 0.0;
    for (int p = 0; p < allPaths.size(); p++) { //2-OPT IS APPLIED FOR EVERY PATH
        int pathNodeNumber = allPaths[p]->nodes.size();
        while (cost < 0.0) {//WHILE COST LESS THAN ZERO THERE'S AN IMPROVEMENT IN ENERGY COST
            cost = 0.0;
                std::list<Node*>::iterator firstNode = allPaths[p]->nodes.begin();
                for (int i = 0; i < pathNodeNumber-1 && cost > 0; i++) { 
                    std::list<Node*>::iterator secondNode = firstNode;
                    secondNode++;
                    for (int j = i + 1; j < pathNodeNumber && cost > 0; j++) {
                        cost = 0.0;
                        if (i == 0) {
                            cost = cost - sqrt((*firstNode)->x * (*firstNode)->x + (*firstNode)->y*(*firstNode)->y) + sqrt((*secondNode)->x*(*secondNode)->x);
                        }
                        else {
                            std::list<Node*>::iterator firstNodePrev = firstNode;
                            firstNodePrev--;
                            cost = cost - travelCosts[(*firstNodePrev)->id][(*firstNode)->id] + travelCosts[(*firstNodePrev)->id][(*secondNode)->id];
                        }
                        if (j == pathNodeNumber-1) {
                            cost = cost - sqrt((*secondNode)->x*(*secondNode)->x + (*secondNode)->y*(*secondNode)->y) + sqrt((*firstNode)->x*(*firstNode)->x + (*firstNode)->y*(*firstNode)->y);
                        }
                        else {
                            std::list<Node*>::iterator secondNodeNext = secondNode;
                            secondNodeNext++;
                            cost = cost - travelCosts[(*secondNode)->id][(*secondNodeNext)->id] + travelCosts[(*firstNode)->id][(*secondNodeNext)->id];
                        }
                        std::list<Node*>::iterator secondNodePrev = secondNode;
                        secondNodePrev--;
                        std::list<Node*>::iterator firstNodeNext = firstNode;
                        firstNodeNext++;
                        cost = cost - travelCosts[(*firstNode)->id][(*firstNodeNext)->id] + travelCosts[(*secondNode)->id][(*firstNodeNext)->id];
                        cost = cost - travelCosts[(*secondNodePrev)->id][(*secondNode)->id] + travelCosts[(*secondNodePrev)->id][(*firstNode)->id];
                        if (cost < 0.0) {
                            std::list<Node*> copy;
                            copy.assign(firstNode, secondNode);
                            copy.reverse();
                            allPaths[i]->nodes.erase(firstNode, secondNode);
                            allPaths[i]->nodes.insert(firstNode, copy.begin(), copy.end());
                            allPaths[i]->energySpent = allPaths[i]->energySpent + cost;
                        //THE PART WHERE THE NEW SCORE IS COMPUTED IS IGNORED AS THE GRAPH IS SYMMETRIC SO THE ONLY CHANGES TO THE ENERGY COST ARE DERIVED FROM THE CHANGES LOCAL TO firstNode AND secondNode. THE CHANGE IS ALREADY CONSIDERED IN THE VARIABLE cost
                        }
                    }
                    firstNode++;
                }
        }
    }
}

//THIS FUNCTION DROPS THE NODE TO A PATH IN NTOP OR A NEW PATH, NOTICE THAT IT'S POSSIBLE FOR THE NEW PATH TO BECOME A TOP PATH
void dropPart(Node* toDrop, std::pair<Path*, std::pair<double, double>>* part) {
    //NOTE: WE ASSUME THAT A mergeNode WAS CALLED BEFORE LAUNCHING THE dropPart
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Path*>::iterator firstNTopPath = tmpAllPaths.begin();
    std::advance(firstNTopPath, M);
    std::vector<Path*> nTopPaths (firstNTopPath, tmpAllPaths.end());
    std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(toDrop, nTopPaths);
    if (result.first.first != -1) {
        movePart(part, result.first.first, result.second.first, result.second.second);
    }
    else {
        Path* pathPtr = new Path(allPaths.size());
        movePart(part, toDrop->id, allPaths.size()-1, 0);
    }
}

//THIS FUNCTION IMPLEMENTS THE REINITIALIZATION I PROCEDURE AS OUTLINED IN LITERATURE (THIS ONE TAKES IN ACCOUNT NODES BASED ON THEIR SCORE)
void reinitializationI(int k) {
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Node*> topNodes;
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> topParts;//THE ERROR IS HERE!!! WE DO NOT USE A POINTER TO THE PART (SECOND MEMBER OF THE EXTERNAL PAIR IN THIS VECTOR) SO THE CHANGES DO NOT APPLY TO THE ORIGINAL PART, INSTEAD THEY APPLY TO A COPY
    for (int i = 0; i < M && i < allPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes[i]->id, &topNodes[i]->parts.pathAssignments[j]));
        }
    }
    std::sort(topParts.begin(), topParts.end(), comparePartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts[i].first < allNodes.size()); i++) {
        dropPart(allNodes[topParts[i].first], topParts[i].second);
    }
    record = calculateRecord(allPaths);
    deviation = (1.0-(prob/100))*record;
}

//THIS FUNCTION IMPLEMENTS THE REINITIALIZATION II PROCEDURE AS OUTLINED IN LITERATURE (THIS ONE TAKES IN ACCOUNT NODES BASED ON THEIR SCORE WEIGHTED ON THE COST)
void reinitializationII(int k){
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Node*> topNodes;
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> topParts;
    for (int i = 0; i < M && i < tmpAllPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes[i]->id, &topNodes[i]->parts.pathAssignments[j]));
        }
    }
    std::sort(topParts.begin(), topParts.end(), compareWeightedPartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts[i].first < allNodes.size()); i++) {
        dropPart(allNodes[topParts[i].first], topParts[i].second);
    }
    record = calculateRecord(allPaths);
    deviation = (1.0-(prob/100))*record;
}

//THIS FUNCTION IMPLEMENTS THE TWO POINT EXCHANGE AS IN LITERATURE
void twoPointExchange() {
    int actionDone = 0;
    int candPath = 0;
    double totScore = 0.0;
    double candScore;
    double cost;
    int topPathPos;
    int candNode;
    int candDestPath;
    int candPositionForTopNode;
    int candPositionForNTopNode;
    std::pair<Path*, std::pair<double, double>>* candPart;
    //CREATE A VECTOR OF THE NODES IN TOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> topNodes;
    //WE CREATE A SORTED VERSION OF allPaths
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    //CREATE A VECTOR OF REFERENCES TO THE NODES IN THE FIRST M PATH (THE TOP ONES)
    for (int i = 0; i < M; i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS OF THE TOP NODES, SAVED IN THE topNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> topParts;
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes[i]->id, &(topNodes[i]->parts.pathAssignments[j])));
            topParts[j].first = topNodes[i]->id;
            topParts[j].second = &(topNodes[i]->parts.pathAssignments[j]);
        }
    }
    //CREATE A VECTOR OF THE NODES IN NTOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> nTopNodes;
    for (int i = M; i < tmpAllPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            nTopNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS IN NTOP PATHS, SAVED IN THE nTopNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> nTopParts;
    for (int i = 0; i < nTopNodes.size(); i++) {
        for (int j = 0; j < nTopNodes[i]->parts.pathAssignments.size(); j++) {
            nTopParts.push_back(std::make_pair(nTopNodes[i]->id, &(nTopNodes[i]->parts.pathAssignments[j])));
            nTopParts[j].first = nTopNodes[i]->id;
            nTopParts[j].second = &(nTopNodes[i]->parts.pathAssignments[j]);
        }
    }
    for (int i = 0; i < topParts.size(); i++) {//FOR EACH PART IN topParts
        actionDone = 0;
        totScore = -DBL_MAX;
        candScore = -DBL_MAX;
        //REMOVE THE PART FROM THE TOP PATH IN ORDER TO GET AN ACCURATE ESTIMATE OF THE COST BY USING fitMove and cheapestInsertion
        Path* topOrigPath = topParts[i].second->first;
        std::list<Node*>::iterator topOrigPositionIt = find(topOrigPath->nodes.begin(), topOrigPath->nodes.end(), allNodes[topParts[i].first]);
        int topOrigPosition = distance(topOrigPath->nodes.begin(), topOrigPositionIt);
        movePart(topParts[i].second, topParts[i].first, -1, 0);
        for (int j = 0; (j < nTopParts.size()) && (actionDone == 0); j++) {//FOR EACH PART IN nTopParts
            if (topParts[i].first != nTopParts[j].first) {//CHECK THAT THE TWO PARTS ARE NOT OF THE SAME NODE
                //CHECK IF THE EXCHANGE IS POSSIBLE, NOTE THAT WE NEED TO CHECK ONLY THE TOP NODE SIDE BECAUSE WE CAN ALWAYS CREATE A NEW NTOP PATH
                std::pair<double, int> toTop = fitMove(allNodes[nTopParts[j].first], topOrigPath, nTopParts[j].second->second.first);
                cost = toTop.first;
                topPathPos = toTop.second;
                if (cost + topOrigPath->energySpent < energy) {//IF THE NTOP PART FITS WE PROCEED TO PERFORM THE TWO movePart
                    //MOVE NTOP PART IN TOP PATH
                    Path* nTopOrigPath = nTopParts[j].second->first;
                    std::list<Node*>::iterator nTopOrigPositionIt = find(nTopOrigPath->nodes.begin(), nTopOrigPath->nodes.end(), allNodes[topParts[i].first]);
                    int nTopOrigPosition = std::min((int) distance(nTopOrigPath->nodes.begin(), nTopOrigPositionIt), 0);
                    movePart(nTopParts[j].second, nTopParts[j].first, topOrigPath->id, topPathPos);
                    //MOVE TOP PART IN NTOP PATH
                    std::vector<Path*> tmpAllPaths = allPaths;
                    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
                    std::vector<Path*>::iterator firstNTopPath = tmpAllPaths.begin();
                    std::advance(firstNTopPath, M);
                    std::vector<Path*> nTopPaths (firstNTopPath, tmpAllPaths.end());
                    std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(allNodes[topParts[i].first], nTopPaths, topParts[i].second);
                    if (result.second.first != -1) {
                        movePart(topParts[i].second, topParts[i].first, result.second.first, result.second.second);
                    }
                    else {
                        Path* newNTopPath = new Path(allPaths.size());
                        movePart(topParts[i].second, topParts[i].first, newNTopPath->id, result.second.second);
                    }
                    totScore = calculateRecord(allPaths);
                    if (totScore > record) {
                        actionDone = 1;
                        nTopParts.erase(std::remove(nTopParts.begin(), nTopParts.end(), nTopParts[j]), nTopParts.end());
                    }
                    else {
                        if (totScore > candScore) {//SAVE BEST MOVE
                            candScore = totScore;
                            candPart = nTopParts[j].second;
                            candNode = nTopParts[j].first;
                            candDestPath = result.second.first;
                            candPositionForTopNode = result.second.second;
                            candPositionForNTopNode = topPathPos;
                        }
                        //REVERT THE TWO movePart
                        movePart(topParts[i].second, topParts[i].first, -1, 0);
                        movePart(nTopParts[j].second, nTopParts[j].first, nTopOrigPath->id, nTopOrigPosition);
                    }
                }
            }
        }
        if (actionDone == 0) {
            //MOVE THE TOP PART BACK INTO ITS PATH
            if (candScore > deviation) {
                //movePart(std::pair<Path*, std::pair<double, double>>* part, int nodeId, int destPath, int position)
                movePart(candPart, candNode, topOrigPath->id, candPositionForNTopNode);
                if (candDestPath == -1) {
                    Path* tmpPath = new Path(allPaths.size());
                    movePart(topParts[i].second, topParts[i].first, tmpPath->id, 0);
                }
                else {
                    movePart(topParts[i].second, topParts[i].first, candDestPath, candPositionForTopNode);
                }
            }
            else {
                movePart(topParts[i].second, topParts[i].first, topOrigPath->id, topOrigPosition);
            }
        }
    }
}

int main(int argc, char** argv) {
//std::cout << "we enter the initializeNodes function" << std::endl;
    initializeNodes(argv[1]);
//std::cout << "we exit the initializeNodes function" << std::endl;
    double p = 5.0;
    int kLim = 5;
    int iLim = 10;
    int stepsNoImpr = 0;
//std::cout << "we enter the initializePaths function" << std::endl;
    initializePaths(allNodes);
//std::cout << "we exit the initializePaths function" << std::endl;
//check();
//std::cout << "we enter the calculateRecord function" << std::endl;
    double record = calculateRecord(allPaths);
//std::cout << "we exit the calculateRecord function" << std::endl;
//check();
    double prevRecord = record;
    deviation = (1.0-(prob/100))*record;
    for(int k = 0; k < kLim; k++) {
        for(int i = 0; i < iLim; i++) {
            int actiondone = 0;
//std::cout << "we enter the twoPointExchange function" << std::endl;
            twoPointExchange();
//std::cout << "we exit the twoPointExchange function" << std::endl;
//check();
//std::cout << "we enter the onePointMove function" << std::endl;
            onePointMove();
//std::cout << "we exit the onePointMove function" << std::endl;
//check();
//std::cout << "we enter the cleanUp function" << std::endl;
            cleanUp();
//std::cout << "we exit the cleanUp function" << std::endl;
//check();
            if (prevRecord == record) {
                break;
            }
        }
        if (record <= prevRecord) {
            stepsNoImpr++;
        }
        else {
            stepsNoImpr = 0;
        }
        if (stepsNoImpr >= 5) {
            break;
        }
//std::cout << "we enter the reinitializationI function" << std::endl;
        reinitializationI(k);
//std::cout << "we exit the reinitializationI function" << std::endl;
//check();
    }
    p = 2.5;
//std::cout << "we enter the reinitializationII function" << std::endl;
    reinitializationII(kLim);
//std::cout << "we exit the reinitializationII function" << std::endl;
//check();
    for(int k=0; k<kLim; k++) {
        for(int i=0; i<iLim; i++) {
            int actiondone = 0;
            twoPointExchange();
            onePointMove();
            cleanUp();
            if (prevRecord == record) break;
        }
        if (record <= prevRecord) {
            stepsNoImpr++;
        }
        else {
            stepsNoImpr = 0;
        }
        if (stepsNoImpr >= 5) break;
        reinitializationI(k);
    }
    //TODO write best solution to file;
check();
}
=======
// this program will reproduce the execution of Chao's algorithm modified to take into account the split visit

#define _USE_MATH_DEFINES

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <list>
#include <algorithm>
#include <float.h>
/*TODO:
fix all the consequent methods: calculateRecord()
*/
int aaa;
double agentSpeed = 1.0;
double sampleRadius = 1.0;
std::vector<std::vector<double> > travelCosts;
class Node;
class Path;
class LocalDist {
    public:
        double stdDeviation;
        
        //THIS FUNCTION RETURNS THE PROBABILITY TO FIND THE TARGET IN THE GIVEN LOCAL DISTRIBUTION MODEL GIVEN THE AMOUNT OF ENERGY/TIME SPENT IN THE LOCAL SEARCH
        double visitProb(double time) {
            double a = sqrt((agentSpeed * (sampleRadius*2) * time) / M_PI);//DISTANCE FROM CENTRAL POINT OF THE DISTRIBUTION
            return erf(a / (stdDeviation * sqrt(2))); //notice that we don't use the (1/2) * (1+erf(x)) form of the CDF because we consider only the distance from the central point
        }
};

class Partitions {
    public:
        std::vector<std::pair<Path*, std::pair<double, double>>> pathAssignments;//THE FIRST FIELD IS THE PATH, THE SECOND IS THE ENERGY SPENT AND THE EXPECTED SCORE
        Node* node;
};

class Node {
    public:
        int id; //[0-M] value equal to the cell in the cost array that refers to this node
        double score;
        double x; //x coordinate
        double y; //y coordinate
        Partitions parts;//how the node visit is split among paths
        LocalDist model;
        
        Node(int nodeId, double nodeX, double nodeY, double nodeScore) {
            id = nodeId;
            x = nodeX;
            y = nodeY;
            score = nodeScore;
            model = LocalDist();
            model.stdDeviation = 10*sampleRadius;
        }
        
        double visitCost() {
            return (model.stdDeviation*model.stdDeviation*M_PI)/(sampleRadius*sampleRadius*M_PI);
        }
};

class Path {
    public:
        std::list<Node*> nodes;
        double totalScore;
        double energySpent;
        int id;
        
        Path(int pathId) {
            id = pathId;
            totalScore = 0;
            energySpent = 0;
        }
};

//VARIABILI GLOBALI
std::vector<Node*> allNodes;
std::vector<Path*> allPaths;
double energy;
double deviation;
double record;
int M;
double prob = 5.0;


bool compareTopVsNTop(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.first->id < M);
}

bool comparePartsByEnergySpent(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.second.first > b.second.first);
}

bool comparePartsByEnergyAvailable(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return ((energy - a.first->energySpent + a.second.first) > (energy - b.first->energySpent + b.second.first));
}

bool comparePartsScore(std::pair<Path*, std::pair<double, double>> a, std::pair<Path*, std::pair<double, double>> b) {
    return (a.second.second > b.second.second);
}

//THIS FUNCTION FIXES THE SCORES OF THE PARTS SO THAT THE SCORES IN TOP PATHS ARE MAXIMIZED AND THE BIGGEST PARTS FAVORED (AS IF THE VISIT OF THOSE WAS EXECUTED FIRST)
void equalizeNode(Node* node) {
    int size = node->parts.pathAssignments.size();
    std::sort(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), compareTopVsNTop);//ORDER PARTS SO THAT TOP PARTS COME FIRST
    int topParts;
    for (int i = 0; i < size; i++) {
        if (node->parts.pathAssignments[i].first->id < M) {
            topParts++;
        }
        else {
            break;
        }
    }
    std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = node->parts.pathAssignments.begin();
    double energyUsed;
    if (topParts > 0) {
        std::advance(it, topParts-1);//THIS ITERATOR POINTS TO THE LAST TOP PATH
        std::sort(node->parts.pathAssignments.begin(), it, comparePartsByEnergySpent);//ORDER TOP PARTS BY ENERGY
    }
    for (int i = 0; i < topParts; i++) {//UPDATE OF THE NTOP PARTS SCORES
        node->parts.pathAssignments[i].second.second = (node->model.visitProb(energyUsed + node->parts.pathAssignments[i].second.first) - node->model.visitProb(energyUsed))*node->score;
        energyUsed = energyUsed + node->parts.pathAssignments[i].second.first;
    }
    it++;//NOW IT POINTS TO THE FIRST NTOP PATH
    std::sort(it, node->parts.pathAssignments.end(), comparePartsByEnergySpent);//ORDER NTOP PARTS BY ENERGY
    for (int i = topParts; i < size; i++) {//UPDATE OF THE NTOP PARTS SCORES
        node->parts.pathAssignments[i].second.second = (node->model.visitProb(energyUsed + node->parts.pathAssignments[i].second.first) - node->model.visitProb(energyUsed))*node->score;
        energyUsed = energyUsed + node->parts.pathAssignments[i].second.first;
    }
}

//TOFIX REFERENCES TO PAIRS
//THIS FUNCTION SPLITS A CERTAIN PARTITION ASSIGNED TO firstPath IN ORDER TO ASSIGN A PART TO secondPath, IT MAXIMIZES THE PART IN secondPath
void splitNode(Node* node, Path* firstPath, Path* secondPath) {
    //TODO: ADD THE CASE WHERE A PART ALREADY EXISTS IN THE TARGET PATH
    std::pair<Path*, std::pair<double, double>> firstPathPart;
    std::pair<Path*, std::pair<double, double>> secondPathPart;
    for (int i = 0; i < node->parts.pathAssignments.size(); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == firstPath) {
            firstPathPart = node->parts.pathAssignments[i];
        }
        if (node->parts.pathAssignments[i].first == secondPath) {
            secondPathPart = node->parts.pathAssignments[i];
        }
    }
    if (secondPathPart.first == NULL) {
        secondPathPart = std::make_pair(secondPath, std::make_pair(0.0, 0.0));
    }
    firstPathPart.second.first = firstPathPart.second.first - std::max(0.0, firstPathPart.second.first - (energy - secondPath->energySpent));
    firstPathPart.first->energySpent = firstPathPart.first->energySpent - std::max(0.0, firstPathPart.second.first - (energy - secondPath->energySpent));
    secondPathPart.second.first = std::min((energy - secondPath->energySpent), firstPathPart.second.first);
    if (firstPathPart.second.first == 0.0) {
        node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), firstPathPart), node->parts.pathAssignments.end());
    }
    equalizeNode(node);       
}

//THIS FUNCTION REVERTS A SPLIT NODE
void revertSplitNode(Node* node, Path* firstPath, Path* secondPath) {
    std::pair<Path*, std::pair<double, double>> firstPathPart;
    std::pair<Path*, std::pair<double, double>> secondPathPart;
    int firstPathFlag = 0;
    int secondPathFlag = 0;
    for (int i = 0; (i < node->parts.pathAssignments.size()) && ((firstPathFlag == 1) && (secondPathFlag == 1)); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == firstPath) {
            firstPathPart = node->parts.pathAssignments[i];
            firstPathFlag = 1;
        }
    }
    for (int i = 0; i < node->parts.pathAssignments.size(); i++) {//FIND THE PART OF firstPath
        if (node->parts.pathAssignments[i].first == secondPath) {
            secondPathPart = node->parts.pathAssignments[i];
            secondPathFlag = 1;
        }
    }
    firstPathPart.second.first = firstPathPart.second.first + secondPathPart.second.first;
    firstPathPart.second.second = firstPathPart.second.second + secondPathPart.second.second;
    node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), secondPathPart), node->parts.pathAssignments.end());
    equalizeNode(node);
}

//this function takes a pointer to a node and returns the cost to reach it and visit it from the start point
double travelCost(Node* a) {
    double x = a->x;
    double y = a->y;
    return sqrt(x*x+y*y)/* + visitCost(a)*/;
}

//this function takes a pointer to a start node "a" and a pointer to a destination node "b" returns the cost to reach it and visit it from the start node
double travelCost(Node* a, Node* b) {
    double x = a->x - b->x;
    double y = a->y - b->y;
    return sqrt(x * x + y * y)/* + visitCost(b)*/;
}

//this function returns the energy cost to visit the node
double visitCost(Node* b) {
    return (b->model.stdDeviation*2)/(sampleRadius*2);
}

//this function takes a vector of pointers to nodes and returns a vector matrix representing the travel+visit costs ordered by node ID number
std::vector<std::vector<double> > updateCosts(std::vector<Node *> allNodes) {
    std::vector<Node *>::iterator it = allNodes.begin();
    //resize the matrix accordingly to the number of nodes
    std::vector<std::vector<double> > distance (allNodes.size());
    for (int i = 0; i < allNodes.size(); i++) {
        distance[i].resize(allNodes.size());
    }
    //update the matrix
    for (int i = 0; i < allNodes.size(); i++) {
        for (int j = 0; j < allNodes.size(); j++) {
            distance[i][j] = travelCost(allNodes[i], allNodes[j]);
            distance[i][j] = distance[i][j];
        } 
    }
    //return statement
    return distance;
}

//this function initializes the nodes by reading their position and scores from a file, it initializes a vector allNodes containing a reference to all the nodes considered.
void initializeNodes(std::string filename) {
    //open file
    std::ifstream myReadFile;
    myReadFile.open(filename);
    //create variables
    std::string line;
    std::string tmp;
    int nodeId = 0;
    double x;
    double y;
    double score;
    Node* newNode;
    if (myReadFile.is_open()) {
        //read line by line
        std::getline(myReadFile, line);
        //divide line in single words
        std::istringstream is( line );
        is >> energy;
        is >> M;
        //first line is the energy available to each agent
        while (!myReadFile.eof()) {
            //the subsequent lines are of the form [x, y, score] where x is the position on the x-axis, y is the position on the y-axis and score is the score associated with the target
            std::getline(myReadFile, line);
            if (myReadFile.eof()) {
                break;
            }
            std::istringstream is( line );
            is >> x;
            is >> y;
            is >> score;
            //create a new node and if its reacheable by the agents add its reference to the allNodes vector, notice that the IDs of nodes are sequential.
            newNode = new Node(nodeId, x, y, score);
            if (travelCost(newNode) < energy) {
                allNodes.push_back(newNode);
                nodeId++;
            }
        }
        //update the travel cost matrix
        travelCosts = updateCosts(allNodes);
    }
    myReadFile.close();
}

//comparison function between nodes to discern which is further from the origin
bool compOrigDistance(Node* a, Node* b) {
    return (sqrt(a->x * a->x + a->y * a->y)) > (sqrt(b->x * b->x + b->y * b->y));
}

//comparison function between paths to discern which has the greater complexive score
bool compPathScore(Path* a, Path* b) {
    return (a->totalScore > b->totalScore);
}


void check() {
    std::cout << "Result:" << std::endl;
    for (int i = 0; i < allPaths.size(); i++) {
        std::cout << "path " << allPaths[i]->id << ": " << std::endl;
        for (std::list<Node*>::iterator it = allPaths[i]->nodes.begin(); it != allPaths[i]->nodes.end(); it++) {
            std::cout << + (*it)->id;
            for (int j = 0; j < (*it)->parts.pathAssignments.size(); j++) {
                if ((*it)->parts.pathAssignments[j].first == allPaths[i]) {
                    std::cout << "(" << (*it)->parts.pathAssignments[j].second.first << "-" << (*it)->parts.pathAssignments[j].second.second << ")";
                    std::cout << " ";
                    break;
                }
                std::cout << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cin >> aaa;
}

//DONE
//TODO: interrogarsi se non sia possibile avere più posizioni ottime di cui una sia preferibile a causa della possibilità di replanning.
//this function takes a pointer to a node, a pointer to a path and returns a pair composed by the minimal cost of inserting the cost into the path and the optimal insertion position.
std::pair<double, int> fitMove(Node* a, Path* path, double visitCost) {
    double minCost = DBL_MAX;
    double cost = DBL_MAX;
    int position = -1;
    int actual = 0;
    int admissible = 0; //FLAG VARIABLE, 1 == AN ADMISSIBLE SOLUTION WAS ALREADY FOUND, 0 == NO ADMISSIBLE SOLUTION WAS FOUND
    if (path->nodes.size() == 0) {
        cost = travelCost(a)*2 + visitCost;
        if (((cost + path->energySpent) < energy) && ((cost + path->energySpent) < energy)) {//NOTE: cost < minCost IS ALWAYS TRUE IN THIS CASE, SO IT'S OMITTED
            //NOTE THERE'S ONLY ONE POSSIBLE POSITION WHEN path->nodes.size() == 0 SO THE POSITION ASSIGNMENT IS OMITTED
            minCost = cost;
            position = actual;
            admissible = 1;
        }
        return std::make_pair(minCost, 0);
    }
    else {
        //POSITION 0
        std::list<Node*>::iterator it = path->nodes.begin();
        std::list<Node*>::iterator itNext = path->nodes.begin();
        itNext++;
        cost = travelCost(a) + travelCosts[a->id][(*it)->id] - travelCost(*it) + visitCost;
        if (cost + path->energySpent < energy && cost + path->energySpent < energy) {
            minCost = cost;
            position = actual;
            admissible = 1;
        }
        while (itNext != path->nodes.end()) {//POSITION 1 TO LAST
            actual++;
            cost = travelCosts[(*it)->id][a->id] + travelCosts[a->id][(*itNext)->id] - travelCosts[(*it)->id][(*itNext)->id] + visitCost;            
            if (cost < minCost) {
                if (admissible == 0) {//NOTE: IN MOST CASES THE CHECK ON THE ENERGY IS NOT NECESSARY SO WE ADDED THE admissible FLAG TO EXPLOIT LAZY EVALUATION THE NESTED IF IS NECESSARY TO SKIP ASSIGNMENTS IN CASE THE SOLUTION IS UNACCEPTABLE
                    if (cost + path->energySpent < energy) {
                        minCost = cost;
                        position = actual;
                        admissible = 1;
                    }
                }
                else {
                    minCost = cost;
                    position = actual;
                }
            }
            it++;
            itNext++;
        }
        //LAST POSITION
        cost = travelCosts[(*it)->id][a->id] + travelCost(a) - travelCost(*it);
        if (cost < minCost) {
            if (admissible == 0) {//NOTE: IN MOST CASES THE CHECK ON THE ENERGY IS NOT NECESSARY SO WE ADDED THE admissible FLAG CHECK TO SPARE SOEM CYCLES
                if (cost + path->energySpent < energy) {
                    minCost = cost;
                    position = actual;
                    //admissible = 1; //NOTE: IT'S NONSENSE WASTING TIME TO SET THE FLAG WHEN WE ALREADY REACHED THE END OF THE PATH
                }
            }
            else {
                minCost = cost;
                position = actual;
            }
        }
        return std::make_pair(minCost, position);
    }
}

//this function takes a vector of nodes, a vector of paths and returns the cheapest insertion of a node from the first vector to a path in the second in the form of pair of pairs. The first pair is the chosen node's id and the cost. The second one is the path and the position where to insert.
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = -1;
    int next = -1;
    double cost;
    int admissible;
    int samePath;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//CHECK ON EACH POSITION
            admissible = 0;
            samePath = 0;
            for (int l = 0; l < nodes[i]->parts.pathAssignments.size(); l++) {//FOR EACH PART BELONGING TO THE NODE 
                if ((nodes[i]->parts.pathAssignments[l].first = paths[j]) && (nodes[i]->parts.pathAssignments[l].first != NULL)) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {
                std::pair<double, int> result = fitMove(nodes[i], paths[j], visitCost(nodes[i]));//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if (cost < minCost && paths[j]->energySpent + cost < energy) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths[j]).id;
                        minNode = nodes[i]->id;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths[j]->id;
                        minNode = nodes[i]->id;
                    }
                }
            }
        }
    }
    return std::make_pair(std::make_pair(minNode, minCost), std::make_pair(minPath, position));
}

//this function takes a vector of nodes, a vector of paths and returns the cheapest insertion of a node from the first vector to a path in the second in the form of pair of pairs. The first pair is the chosen node's id and the cost. The second one is the path and the position where to insert.
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths, std::pair<Path*, std::pair<double, double>>* part) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = -1;
    int next = -1;
    int samePath = 0;;
    double cost;
    int admissible;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//CHECK ON EACH POSITION
            admissible = 0;
            samePath = 0;
            for (std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = nodes[i]->parts.pathAssignments.begin(); it != nodes[i]->parts.pathAssignments.end(); it++) {
                if ((*it).first = paths[j]) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {
                std::pair<double, int> result = fitMove(nodes[i], paths[j], part->second.first);//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if (cost < minCost && paths[j]->energySpent + cost < energy) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths[j]).id;
                        minNode = i;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths[j]->id;
                        minNode = nodes[i]->id;
                    }
                }
            }
        }
    }
    return std::make_pair(std::make_pair(minNode, minCost), std::make_pair(minPath, position));
}

//THIS FUNCTION TAKES A REFERENCE TO THE PART MOVED, AN INTEGER REPRESENTING THE ID OF THE NODE, AN INT REPRESENTING THE ID OF THE DESTINATION PATH AND AN INT WHICH IS THE POSITION IN THE NEW PATH. THEN IT MOVES THE NODE WITH SUCH AN ID IN THE POSITION IN THE NEW PATH AND UPDATES THE path FIELD IN THE NODE
void movePart(std::pair<Path*, std::pair<double, double>>* part, int nodeId, int destPath, int position) {
    if (part->first != NULL) {//IF THE PREVIOUS PATH IS NOT NULL REMOVE PART FROM THE PREVIOUS PATH
        Path* originPath = part->first;
        part->first = NULL;
        std::list<Node*>::iterator it = find(originPath->nodes.begin(), originPath->nodes.end(), allNodes[nodeId]);
        originPath->totalScore = originPath->totalScore - part->second.second;
            std::list<Node*>::iterator itNext = it;
            itNext++;
        if (it != originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCosts[(*itPrev)->id][nodeId] - travelCosts[nodeId][(*itNext)->id] + travelCosts[(*itPrev)->id][(*itNext)->id];
        }
        if (it == originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCost(allNodes[nodeId]) - travelCosts[nodeId][(*itNext)->id];
        }
        if (it != originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            originPath->energySpent = originPath->energySpent - travelCosts[(*itPrev)->id][nodeId] - travelCost(allNodes[nodeId]);
        }
        if (it == originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            originPath->energySpent = originPath->energySpent - (travelCost(allNodes[nodeId])*2);
        }
        originPath->nodes.erase(std::remove(originPath->nodes.begin(), originPath->nodes.end(), allNodes[nodeId]), originPath->nodes.end());
    }
    //ADD PART TO THE NEW PATH
    if (destPath >= 0) {//IF THE id OF THE DESTINATION PATH IS POSITIVE THE DESTINATION PATH IS ASSUMED TO EXIST
        std::list<Node*>::iterator it = allPaths[destPath]->nodes.begin();
        std::advance(it, position);
        part->first = allPaths[destPath];
        allPaths[destPath]->nodes.insert(it, allNodes[nodeId]);
        it = allPaths[destPath]->nodes.begin();
        std::advance(it, position);
        allPaths[destPath]->totalScore = allPaths[destPath]->totalScore + part->second.second;
        if (position != 0 && position != allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCosts[(*itPrev)->id][nodeId] + travelCosts[nodeId][(*itNext)->id];
        }
        if (position == 0 && position != allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCost(allNodes[nodeId]) + travelCosts[nodeId][(*itNext)->id];
        }
        if (position != 0 && position == allPaths[destPath]->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCosts[(*itPrev)->id][nodeId] + travelCost(allNodes[nodeId]);
        }
        if (position == 0 && position == allPaths[destPath]->nodes.size()-1) {
            allPaths[destPath]->energySpent = allPaths[destPath]->energySpent + travelCost(allNodes[nodeId])*2;
        }
    }
}

//THIS FUNCTION CONSIDERS THE PARTITIONS OF A NODE AND TRIES TO MERGE THE PARTS IN A NTOP PATH TO THE PARTS IN A TOP PATH
void mergeNode(Node* node) {
    int size = node->parts.pathAssignments.size();
    std::sort(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), compareTopVsNTop);//ORDER PARTS SO THAT TOP PARTS COME FIRST
    int topParts = 0;
    for (int i = 0; i < size; i++) {
        if (node->parts.pathAssignments[i].first->id < M) {
            topParts++;
        }
        else {
            break;
        }
    }
    std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = node->parts.pathAssignments.begin();
    if (topParts > 0) {
        std::advance(it, topParts-1);// std::(it, topParts-1);//THIS ITERATOR POINTS TO THE LAST TOP PART
        std::sort(node->parts.pathAssignments.begin(), it, comparePartsByEnergyAvailable);//NOTE: ORDER TOP PARTS BY SUM OF ENERGY LEFT IN THE PATH AND ENERGY ALREADY ASSIGNED TO THE VISIT
        it++;//NOW IT POINTS TO THE FIRST NTOP PATH
    }
    std::sort(it, node->parts.pathAssignments.end(), comparePartsByEnergyAvailable);//NOTE: ORDER NTOP PARTS BY SUM OF ENERGY LEFT IN THE PATH AND ENERGY ALREADY ASSIGNED TO THE VISIT
    for (int i = 0; i < topParts; i++) {
        for (int j = size - 1; j > (topParts - 1); j--) {
            node->parts.pathAssignments[i].second.first = node->parts.pathAssignments[i].second.first + std::min((energy - node->parts.pathAssignments[i].first->energySpent), node->parts.pathAssignments[j].second.first);
            node->parts.pathAssignments[j].second.first = std::max(0.0, (node->parts.pathAssignments[j].second.first - (energy - node->parts.pathAssignments[i].first->energySpent)));
            if ((energy - node->parts.pathAssignments[i].first->energySpent) == 0) {
                break;
            }
        }
        for (int j = topParts - 1; j > i; j--) {
            node->parts.pathAssignments[i].second.first = node->parts.pathAssignments[i].second.first + std::min((energy - node->parts.pathAssignments[i].first->energySpent), node->parts.pathAssignments[j].second.first);
            node->parts.pathAssignments[j].second.first = std::max(0.0, (node->parts.pathAssignments[j].second.first - (energy - node->parts.pathAssignments[i].first->energySpent)));
            if ((energy - node->parts.pathAssignments[i].first->energySpent) == 0) {
                break;
            }
        }
    }
    for (int i = 0; i < topParts; i++) {
        if (node->parts.pathAssignments[i].second.first == 0.0) {
            movePart(&node->parts.pathAssignments[i], node->id, -1, 0);
            node->parts.pathAssignments[i].first->nodes.erase(std::remove(node->parts.pathAssignments[i].first->nodes.begin(), node->parts.pathAssignments[i].first->nodes.end(), node), node->parts.pathAssignments[i].first->nodes.end());
            node->parts.pathAssignments.erase(std::remove(node->parts.pathAssignments.begin(), node->parts.pathAssignments.end(), node->parts.pathAssignments[i]), node->parts.pathAssignments.end());
        }
    }
    equalizeNode(node);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, std::vector<Path*> paths, std::pair<Path*, std::pair<double, double>>* part) {
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths, part);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, std::vector<Path*> paths) {
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(Node* node, Path* path) {
    std::vector<Path*> paths;
    paths.push_back(path);
    std::vector<Node*> nodes;
    nodes.push_back(node);
    return cheapestInsertion(nodes, paths);
}

//ADAPTER FUNCTION THAT INCAPSULATE THE PARAMETERS IN VECTORS AND INVOKES cheapestInsertion
std::pair< std::pair<int, double>, std::pair<int, int>> cheapestInsertion(std::vector<Node*> nodes, Path* path) {
    std::vector<Path*> paths;
    paths.push_back(path);
    return cheapestInsertion(nodes, paths);
}

//BOOLEAN FUNCTION TO SORT PATHS BY SCORE (GREATER SCORE FIRST)
bool scorePathComp(Path* a, Path* b) {
    return a->totalScore > b->totalScore;
}

//THIS FUNCTION GETS A COPY OF THE topPaths VECTOR AND SUMS THE SCORES OF THE PATHS IN IT
double calculateRecord(std::vector<Path*> paths){
    std::vector<Path*> tmpAllPaths = paths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    double sum = 0.0;
    for (int i = 0; i < M; i++) {
        sum = sum + tmpAllPaths[i]->totalScore;
    }
    return sum;
}

//THIS FUNCTION CREATES THE INITIAL PATHS. IT ADDS THE FURTHEST NODES TO A FEW PATHS AND THEN ADDS THE REMAINING NODES WITH CHEAPEST INSERTION
void initializePaths(std::vector<Node*> nodesToAdd) {
    int nodesNumber = nodesToAdd.size();
    int L = std::min(5, nodesNumber);//NOTE: IN THE ORIGINAL WORK PRESENTED BY CHAO THE NUMBER OF AGENTS IS LOWER THAN 5. IF WE GO ABOVE WE MAY NEED TO CHANGE THIS LINE OF CODE.
    Path* pathPtr;
    //sort the nodes in descending order of distance from origin
    std::sort(nodesToAdd.begin(), nodesToAdd.end(), compOrigDistance);
    for (int i = 0; i < M && i < L; i++) {//FOR EACH INITIAL PATH (THE NUMBER OF INITIAL PATH IS THE NUMBER OF AGENTS) ADD THE FURTHEST NODE FROM THE START
        pathPtr = new Path(i);
        allPaths.push_back(pathPtr);
        nodesToAdd[i]->parts.pathAssignments.push_back(std::make_pair(pathPtr, std::make_pair(visitCost(nodesToAdd[i]), nodesToAdd[i]->score)));
        movePart(&nodesToAdd[i]->parts.pathAssignments[0], nodesToAdd[i]->id, pathPtr->id, pathPtr->nodes.size());
        nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), nodesToAdd[i]), nodesToAdd.end());
    }
    int prevSize;
    int pathNumber = M;
//INFINITE LOOP BELOW, APPARENTLY THE nodesToAdd SIZE DIMINISHES ONLY IN THE FIRST CYCLE
    while (nodesToAdd.size() > 0) {//WHILE THERE ARE NODES WITHOUT PATH WE INSERT ONE OF THOSE IN A PATH, WHEN NECESSARY WE CREATE A NEW ONE
        prevSize = nodesToAdd.size();
        std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(nodesToAdd, allPaths);
        if (result.first.first == -1) {//IF THERE'S NO ADMISSIBLE INSERTION THEN CREATE A NEW PATH AND REPEAT
            pathPtr = new Path(pathNumber);
            pathNumber++;
            allPaths.push_back(pathPtr);
        }
        else {//ELSE ADD THE NODE TO THE PATH, UPDATE THE NODE'S PATH FIELD AND THE PATH'S SCORE AND ENERGY SPENT FIELDS
            allNodes[result.first.first]->parts.pathAssignments.push_back(std::make_pair(allPaths[result.second.first], std::make_pair(visitCost(allNodes[result.first.first]), allNodes[result.first.first]->score)));
            movePart(&allNodes[result.first.first]->parts.pathAssignments[0], result.first.first, result.second.first, result.second.second);
            nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), allNodes[result.first.first]), nodesToAdd.end());
        }
    }
    record = calculateRecord(allPaths);//SOME NON-TRIVIAL WORK WAS DONE ON THE allPaths VECTOR SO WE UPDATE THE record
}

//ARRIVATO QUI
//NOTE:
//1. DOBBIAMO CONSIDERARE IL COSTO DI VISITA SEMPRE SEPARATO DAL COSTO DI VIAGGIO
//2. QUANDO DECIDIAMO SE SPLITTARE UN NODO DOBBIAMO CONSIDERARE SE NE VALE LA PENA, SUGGERISCO DI IMPOSTARE UNA PERCENTUALE PARI A QUELLA DELLA DEVIATION E SE UNA QUALCHE FORMA DI RAPPORTO SCORE/COSTO SCENDE SOTTO QUELLA PERCENTUALE NON ESEGUIRE LO SPLIT
//TODO (MAYBE): ADD A METHOD THAT UPDATES THE SCORE OF PATHS
//TODO: ADD A GLOBAL VECTOR TOP PATHS AND A GLOBAL VECTOR NTOP PATHS SO THAT THEY'RE ALWAYS ACCOUNTED FOR CORRECTLY
//TODO: ADD A isTopPath(Path* path) BOOLEAN FUNCTION THAT CHECKS IF path IS A TOP PATH AND USE IT INSTEAD OF CHECKING THE PATH ID TO KNOW IF THE PATH IS A TOP PATH
//TODO: FIX THE WHOLE allPaths[] ORDER SITUATION BY LETTING THE calculateRecord() CREATE A COPY OF THE VECTOR TO ORDER AND CALCULATE THE RECORD ON THAT. OR BY ADDING A topPaths and a nTopPaths VECTOR.
//TODO: AGGIUNGERE UN CASO A splitNode NEL CASO IN CUI IL NODO DI CUI FACCIAMO LO SPLIT HA GIA' UNA PARTE NEL DESTINATION PATH

//THIS FUNCTION IMPLEMENTS THE ONE-POINT MOVEMENT AS OUTLINED IN THE LITERATURE
void onePointMove() {
    int bestMovePosition;
    int bestMovePath = -1;
    std::pair<Path*, std::pair<double, double>>* bestMovePart;
    std::pair<Path*, std::pair<double, double>>* candPart;
    double bestMove = 0;
    double cost = 0;
    int actionDone;
    int prevPath;
    int prevPosition;
    double candRecord;
    double totScore; 
    for (int i = 0; i < allNodes.size(); i++) {//FOR EACH NODE
        bestMovePath = -1;
        bestMovePosition = -1;
        candRecord = -DBL_MAX;
        actionDone = 0;
        mergeNode(allNodes[i]);
        for (int l = 0; (l < allNodes[i]->parts.pathAssignments.size()) && (actionDone == 0); l++) {//FOR EACH PART
            candPart = &allNodes[i]->parts.pathAssignments[l];
            prevPath = candPart->first->id;
            prevPosition = 0;
            for (std::list<Node*>::iterator tmpIt = allPaths[prevPath]->nodes.begin(); tmpIt != allPaths[prevPath]->nodes.end(); tmpIt++) {
                if (*(tmpIt) == allNodes[i]) {//IF WE'VE FOUND THE POSITION WE WXIT THE FOR LOOP
                    break;
                }
                prevPosition++;
            }
            for (int j = 0; (j < allPaths.size()) && (actionDone == 0); j++) {//WE TRY TO MOVE IT IN ANOTHER PATH
                Path* candPath = allPaths[j];
                if (candPath->id != prevPath) {
                    std::pair<double, int> result = fitMove(allNodes[i], candPath, candPart->second.first);
                    cost = result.first;
                    if (cost + candPath->energySpent < energy) {
                        movePart(candPart, i, j, result.second);
                        totScore = calculateRecord(allPaths);
                        if (totScore > record) {
                            record = totScore;
                            deviation = (1.0-(prob/100))*record;
                            actionDone = 1;
                        }
                        else {
                        movePart(candPart, i, prevPath, prevPosition);//SEGFAULT HERE, THE LAST PARAMETER IS WRONG
                            if (totScore > deviation && totScore > candRecord) {//IF THE SOLUTION IS NOT A STRICT IMPROVEMENT WE SAVE THE BEST SOLUTION THAT DOES NOT WORSEN THE SOLUTION BY MORE THAN THE DEVIATION
                                bestMovePath = j;
                                bestMovePosition = result.second;
                                candRecord = totScore;
                                bestMovePart = &allNodes[i]->parts.pathAssignments[l];
                            }
                        }
                    }
                    else {
                        //SAVE PATH AND NODE STATUS
                        Node tmpNode = *(allNodes[i]);
                        Path tmpPath1 = *(allPaths[j]);
                        Path tmpPath2 = *(allPaths[prevPath]);
                        //PERFORM SPLIT NODE
                        if ((result.first / candPart->second.second) > (prob/100)) {
                            splitNode(allNodes[i], allPaths[prevPath], allPaths[j]);
                            totScore = calculateRecord(allPaths);
                        }
                        if (totScore > record) {//CONFIRM SPLIT NODE
                            record = totScore;
                            deviation = (1.0-(prob/100))*record;
                            actionDone = 1;
                        }
                        else {//REVERT SPLIT NODE
                            *(allNodes[i]) = tmpNode;
                            *(allPaths[j]) = tmpPath1;
                            *(allPaths[prevPath]) = tmpPath2;
                        }
                    }
                }
            }
        }
        if (bestMovePath != -1 && actionDone == 0) {
            movePart(bestMovePart, i, bestMovePath, bestMovePosition);
            record = calculateRecord(allPaths);
            deviation = (1.0-(prob/100))*record;
        }
    }
}

//COMPARISON FUNCTION BETWEEN NODES BASED ON SCORE
bool compareNodeScore(Node* a, Node* b) {
    if (a->score < b->score) {
        return false;
    }
    else {
        return true;
    }
}

//COMPARISON FUNCTION BETWEEN PARTS BASED ON SCORE
/*
bool comparePartScore(std::pair<Path*, std::pair<double, double>>* a, std::pair<Path*, std::pair<double, double>>* b) {
    if (a->second.second < b->second.second) {
        return false;
    }
    else {
        return true;
    }
}*/

bool comparePartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>> a, std::pair<int, std::pair<Path*, std::pair<double, double>>> b) {
    if (a.second.second.second < b.second.second.second) {
        return false;
    }
    else {
        return true;
    }
}

//COMPARISON FUNCTION BETWEEN NODES BASED ON SCORE WEIGHTED ON THE LOCAL COST OF THAT NODE
/*bool compareWeightedNodeScore(Node* a, Node* b) {
    std::list<Node*>::iterator itA =  find(a->path->nodes.begin(), a->path->nodes.end(), a);
    std::list<Node*>::iterator aNext = itA;
    std::list<Node*>::iterator aPrev = itA;
    std::list<Node*>::iterator itB =  find(b->path->nodes.begin(), b->path->nodes.end(), b);
    std::list<Node*>::iterator bNext = itA;
    std::list<Node*>::iterator bPrev = itA;
    aNext++;
    aPrev--;
    bNext++;
    bPrev--;
    double costA = travelCosts[(*aPrev)->id][a->id] + travelCosts[a->id][(*aNext)->id] - travelCosts[(*aPrev)->id][(*aNext)->id];
    double costB = travelCosts[(*bPrev)->id][a->id] + travelCosts[a->id][(*bNext)->id] - travelCosts[(*bPrev)->id][(*bNext)->id];
    return (a->score/costA < b->score/costB);
}*/


//COMPARISON FUNCTION BETWEEN PARTS BASED ON SCORE WEIGHTED ON THE LOCAL COST OF THAT PART
bool compareWeightedPartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>> a, std::pair<int, std::pair<Path*, std::pair<double, double>>> b) {
    std::list<Node*>::iterator itA =  find(a.second.first->nodes.begin(), a.second.first->nodes.end(), allNodes[a.first]);
    std::list<Node*>::iterator itALast = a.second.first->nodes.end();
    itALast--;
    std::list<Node*>::iterator itB =  find(b.second.first->nodes.begin(), b.second.first->nodes.end(), allNodes[b.first]);
    std::list<Node*>::iterator itBLast = b.second.first->nodes.end();
    itBLast--;
    double costA = 0.0;
    if ((itA != a.second.first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        aNext++;
        double costA = travelCosts[(*aPrev)->id][(*itA)->id] + travelCosts[(*itA)->id][(*aNext)->id] - travelCosts[(*aPrev)->id][(*aNext)->id];
    }
    if ((itA != a.second.first->nodes.begin()) && (itA == itALast)) {
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        double costA = travelCosts[(*aPrev)->id][(*itA)->id] + travelCost(allNodes[(*itA)->id]) - travelCost(allNodes[(*aPrev)->id]);
    }
    if ((itA == a.second.first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        aNext++;
        double costA = travelCost(allNodes[(*itA)->id]) + travelCosts[(*itA)->id][(*aNext)->id] - travelCost(allNodes[(*aNext)->id]);
    }
    if ((itA == a.second.first->nodes.begin()) && (itA == itALast)) {
        double costA = travelCost(allNodes[(*itA)->id]) + travelCost(allNodes[(*itA)->id]);
    }
    double costB = 0.0;
    if ((itB != b.second.first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        bNext++;
        double costB = travelCosts[(*bPrev)->id][(*itB)->id] + travelCosts[(*itB)->id][(*bNext)->id] - travelCosts[(*bPrev)->id][(*bNext)->id];
    }
    if ((itB != b.second.first->nodes.begin()) && (itB == itBLast)) {
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        double costB = travelCosts[(*bPrev)->id][(*itB)->id] + travelCost(allNodes[(*itB)->id]) - travelCost(allNodes[(*bPrev)->id]);
    }
    if ((itB == b.second.first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        bNext++;
        double costB = travelCost(allNodes[(*itB)->id]) + travelCosts[(*itB)->id][(*bNext)->id] - travelCost(allNodes[(*bNext)->id]);
    }
    if ((itB == b.second.first->nodes.begin()) && (itB == itBLast)) {
        double costB = travelCost(allNodes[(*itB)->id]) + travelCost(allNodes[(*itB)->id]);
    }
    return (a.second.second.second/costA < b.second.second.second/costB);
}

//THIS FUNCTION IMPLEMENTS THE CLEAN UP ROUTINE
void cleanUp() {
    double record = 0.0;
    int actionDone = 0;
    double cost = 0.0;
    for (int p = 0; p < allPaths.size(); p++) { //2-OPT IS APPLIED FOR EVERY PATH
        int pathNodeNumber = allPaths[p]->nodes.size();
        while (cost < 0.0) {//WHILE COST LESS THAN ZERO THERE'S AN IMPROVEMENT IN ENERGY COST
            cost = 0.0;
                std::list<Node*>::iterator firstNode = allPaths[p]->nodes.begin();
                for (int i = 0; i < pathNodeNumber-1 && cost > 0; i++) { 
                    std::list<Node*>::iterator secondNode = firstNode;
                    secondNode++;
                    for (int j = i + 1; j < pathNodeNumber && cost > 0; j++) {
                        cost = 0.0;
                        if (i == 0) {
                            cost = cost - sqrt((*firstNode)->x * (*firstNode)->x + (*firstNode)->y*(*firstNode)->y) + sqrt((*secondNode)->x*(*secondNode)->x);
                        }
                        else {
                            std::list<Node*>::iterator firstNodePrev = firstNode;
                            firstNodePrev--;
                            cost = cost - travelCosts[(*firstNodePrev)->id][(*firstNode)->id] + travelCosts[(*firstNodePrev)->id][(*secondNode)->id];
                        }
                        if (j == pathNodeNumber-1) {
                            cost = cost - sqrt((*secondNode)->x*(*secondNode)->x + (*secondNode)->y*(*secondNode)->y) + sqrt((*firstNode)->x*(*firstNode)->x + (*firstNode)->y*(*firstNode)->y);
                        }
                        else {
                            std::list<Node*>::iterator secondNodeNext = secondNode;
                            secondNodeNext++;
                            cost = cost - travelCosts[(*secondNode)->id][(*secondNodeNext)->id] + travelCosts[(*firstNode)->id][(*secondNodeNext)->id];
                        }
                        std::list<Node*>::iterator secondNodePrev = secondNode;
                        secondNodePrev--;
                        std::list<Node*>::iterator firstNodeNext = firstNode;
                        firstNodeNext++;
                        cost = cost - travelCosts[(*firstNode)->id][(*firstNodeNext)->id] + travelCosts[(*secondNode)->id][(*firstNodeNext)->id];
                        cost = cost - travelCosts[(*secondNodePrev)->id][(*secondNode)->id] + travelCosts[(*secondNodePrev)->id][(*firstNode)->id];
                        if (cost < 0.0) {
                            std::list<Node*> copy;
                            copy.assign(firstNode, secondNode);
                            copy.reverse();
                            allPaths[i]->nodes.erase(firstNode, secondNode);
                            allPaths[i]->nodes.insert(firstNode, copy.begin(), copy.end());
                            allPaths[i]->energySpent = allPaths[i]->energySpent + cost;
                        //THE PART WHERE THE NEW SCORE IS COMPUTED IS IGNORED AS THE GRAPH IS SYMMETRIC SO THE ONLY CHANGES TO THE ENERGY COST ARE DERIVED FROM THE CHANGES LOCAL TO firstNode AND secondNode. THE CHANGE IS ALREADY CONSIDERED IN THE VARIABLE cost
                        }
                    }
                    firstNode++;
                }
        }
    }
}

//THIS FUNCTION DROPS THE NODE TO A PATH IN NTOP OR A NEW PATH, NOTICE THAT IT'S POSSIBLE FOR THE NEW PATH TO BECOME A TOP PATH
void dropPart(Node* toDrop, std::pair<Path*, std::pair<double, double>>* part) {
    //NOTE: WE ASSUME THAT A mergeNode WAS CALLED BEFORE LAUNCHING THE dropPart
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Path*>::iterator firstNTopPath = tmpAllPaths.begin();
    std::advance(firstNTopPath, M);
    std::vector<Path*> nTopPaths (firstNTopPath, tmpAllPaths.end());
    std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(toDrop, nTopPaths);
    if (result.first.first != -1) {
        movePart(part, result.first.first, result.second.first, result.second.second);
    }
    else {
        Path* pathPtr = new Path(allPaths.size());
        allPaths.push_back(pathPtr);
        movePart(part, toDrop->id, allPaths.size()-1, 0);
    }
}

//THIS FUNCTION IMPLEMENTS THE REINITIALIZATION I PROCEDURE AS OUTLINED IN LITERATURE (THIS ONE TAKES IN ACCOUNT NODES BASED ON THEIR SCORE)
void reinitializationI(int k) {
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Node*> topNodes;
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>>> topParts;
    for (int i = 0; i < M && i < allPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(i, topNodes[i]->parts.pathAssignments[j]));
        }
    }
    std::sort(topParts.begin(), topParts.end(), comparePartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts[i].first < allNodes.size()); i++) {
        dropPart(allNodes[topParts[i].first], &topParts[i].second);
    }
    record = calculateRecord(allPaths);
    deviation = (1.0-(prob/100))*record;
}

//THIS FUNCTION IMPLEMENTS THE REINITIALIZATION II PROCEDURE AS OUTLINED IN LITERATURE (THIS ONE TAKES IN ACCOUNT NODES BASED ON THEIR SCORE WEIGHTED ON THE COST)
void reinitializationII(int k){
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    std::vector<Node*> topNodes;
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>>> topParts;
    for (int i = 0; i < M && i < tmpAllPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes[i]->id, topNodes[i]->parts.pathAssignments[j]));
        }
    }
    std::sort(topParts.begin(), topParts.end(), compareWeightedPartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts[i].first < allNodes.size()); i++) {
        dropPart(allNodes[topParts[i].first], &topParts[i].second);
    }
    record = calculateRecord(allPaths);
    deviation = (1.0-(prob/100))*record;
}

//ARRIVATO QUI
//NOTE:
//1. DOBBIAMO CONSIDERARE IL COSTO DI VISITA SEMPRE SEPARATO DAL COSTO DI VIAGGIO
//2. QUANDO DECIDIAMO SE SPLITTARE UN NODO DOBBIAMO CONSIDERARE SE NE VALE LA PENA, SUGGERISCO DI IMPOSTARE UNA PERCENTUALE PARI A QUELLA DELLA DEVIATION E SE UNA QUALCHE FORMA DI RAPPORTO SCORE/COSTO SCENDE SOTTO QUELLA PERCENTUALE NON ESEGUIRE LO SPLIT
//TODO (MAYBE): ADD A METHOD THAT UPDATES THE SCORE OF PATHS
//TODO: ADD A GLOBAL VECTOR TOP PATHS AND A GLOBAL VECTOR NTOP PATHS SO THAT THEY'RE ALWAYS ACCOUNTED FOR CORRECTLY
//TODO: ADD A isTopPath(Path* path) BOOLEAN FUNCTION THAT CHECKS IF path IS A TOP PATH AND USE IT INSTEAD OF CHECKING THE PATH ID TO KNOW IF THE PATH IS A TOP PATH
//TODO: FIX THE WHOLE allPaths[] ORDER SITUATION BY LETTING THE calculateRecord() CREATE A COPY OF THE VECTOR TO ORDER AND CALCULATE THE RECORD ON THAT. OR BY ADDING A topPaths and a nTopPaths VECTOR.

//THIS FUNCTION IMPLEMENTS THE TWO POINT EXCHANGE AS IN LITERATURE
void twoPointExchange() {
    int actionDone = 0;
    int candPath = 0;
    double totScore = 0.0;
    double candScore;
    double cost;
    int topPathPos;
    int candNode;
    int candDestPath;
    int candPositionForTopNode;
    int candPositionForNTopNode;
    std::pair<Path*, std::pair<double, double>>* candPart;
    //CREATE A VECTOR OF THE NODES IN TOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> topNodes;
    //WE CREATE A SORTED VERSION OF allPaths
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    //CREATE A VECTOR OF REFERENCES TO THE NODES IN THE FIRST M PATH (THE TOP ONES)
    for (int i = 0; i < M; i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS OF THE TOP NODES, SAVED IN THE topNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> topParts;
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes[i]->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes[i]->id, &(topNodes[i]->parts.pathAssignments[j])));
        }
    }
    //CREATE A VECTOR OF THE NODES IN NTOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> nTopNodes;
    for (int i = M; i < tmpAllPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths[i]->nodes.begin(); it != tmpAllPaths[i]->nodes.end(); it++) {
            nTopNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS IN NTOP PATHS, SAVED IN THE nTopNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> nTopParts;
    for (int i = 0; i < nTopNodes.size(); i++) {
        for (int j = 0; j < nTopNodes[i]->parts.pathAssignments.size(); j++) {
            nTopParts.push_back(std::make_pair(nTopNodes[i]->id, &(nTopNodes[i]->parts.pathAssignments[j])));
        }
    }
    for (int i = 0; i < topParts.size(); i++) {//FOR EACH PART IN topParts
        actionDone = 0;
        totScore = -DBL_MAX;
        candScore = -DBL_MAX;
        //REMOVE THE PART FROM THE TOP PATH IN ORDER TO GET AN ACCURATE ESTIMATE OF THE COST BY USING fitMove and cheapestInsertion
        Path* topOrigPath = topParts[i].second->first;
        std::list<Node*>::iterator topOrigPositionIt = find(topOrigPath->nodes.begin(), topOrigPath->nodes.end(), allNodes[topParts[i].first]);
        int topOrigPosition = distance(topOrigPath->nodes.begin(), topOrigPositionIt);
        movePart(topParts[i].second, topParts[i].first, -1, 0);
        for (int j = 0; (j < nTopParts.size()) && (actionDone == 0); j++) {//FOR EACH PART IN nTopParts
            if (topParts[i].first != nTopParts[j].first) {//CHECK THAT THE TWO PARTS ARE NOT OF THE SAME NODE
                //CHECK IF THE EXCHANGE IS POSSIBLE, NOTE THAT WE NEED TO CHECK ONLY THE TOP NODE SIDE BECAUSE WE CAN ALWAYS CREATE A NEW NTOP PATH
                std::pair<double, int> toTop = fitMove(allNodes[nTopParts[j].first], topOrigPath, nTopParts[j].second->second.first);
                cost = toTop.first;
                topPathPos = toTop.second;
                if (cost + topOrigPath->energySpent < energy) {//IF THE NTOP PART FITS WE PROCEED TO PERFORM THE TWO movePart
                    //MOVE NTOP PART IN TOP PATH
                    Path* nTopOrigPath = nTopParts[j].second->first;
                    std::list<Node*>::iterator nTopOrigPositionIt = find(nTopOrigPath->nodes.begin(), nTopOrigPath->nodes.end(), allNodes[topParts[i].first]);
                    int nTopOrigPosition = std::min((int) distance(nTopOrigPath->nodes.begin(), nTopOrigPositionIt), 0);
                    movePart(nTopParts[j].second, nTopParts[j].first, topOrigPath->id, topPathPos);
                    //MOVE TOP PART IN NTOP PATH
                    std::vector<Path*> tmpAllPaths = allPaths;
                    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
                    std::vector<Path*>::iterator firstNTopPath = tmpAllPaths.begin();
                    std::advance(firstNTopPath, M);
                    std::vector<Path*> nTopPaths (firstNTopPath, tmpAllPaths.end());
                    std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(allNodes[topParts[i].first], nTopPaths, topParts[i].second);
                    if (result.second.first != -1) {
                        movePart(topParts[i].second, topParts[i].first, result.second.first, result.second.second);
                    }
                    else {
                        Path* newNTopPath = new Path(allPaths.size());
                        allPaths.push_back(newNTopPath);
                        movePart(topParts[i].second, topParts[i].first, newNTopPath->id, result.second.second);
                    }
                    totScore = calculateRecord(allPaths);
                    if (totScore > record) {
                        actionDone = 1;
                        nTopParts.erase(std::remove(nTopParts.begin(), nTopParts.end(), nTopParts[j]), nTopParts.end());
                    }
                    else {
                        if (totScore > candScore) {//SAVE BEST MOVE
                            candScore = totScore;
                            candPart = nTopParts[j].second;
                            candNode = nTopParts[j].first;
                            candDestPath = result.second.first;
                            candPositionForTopNode = result.second.second;
                            candPositionForNTopNode = topPathPos;
                        }
                        //REVERT THE TWO movePart
                        movePart(topParts[i].second, topParts[i].first, -1, 0);
                        movePart(nTopParts[j].second, nTopParts[i].first, nTopOrigPath->id, nTopOrigPosition);
                    }
                }
            }
        }
        if (actionDone == 0) {
            //MOVE THE TOP PART BACK INTO ITS PATH
            if (candScore > deviation) {
                //movePart(std::pair<Path*, std::pair<double, double>>* part, int nodeId, int destPath, int position)
                movePart(candPart, candNode, topOrigPath->id, candPositionForNTopNode);
                if (candDestPath == -1) {
                    Path* tmpPath = new Path(allPaths.size());
                    allPaths.push_back(tmpPath);
                    movePart(topParts[i].second, topParts[i].first, tmpPath->id, 0);
                }
                else {
                    movePart(topParts[i].second, topParts[i].first, candDestPath, candPositionForTopNode);
                }
            }
            else {
                movePart(topParts[i].second, topParts[i].first, topOrigPath->id, topOrigPosition);
            }
        }
    }
}

int main(int argc, char** argv) {
std::cout << "we enter the initializeNodes function" << std::endl;
    initializeNodes(argv[1]);
std::cout << "we exit the initializeNodes function" << std::endl;
    double p = 5.0;
    int kLim = 5;
    int iLim = 10;
    int stepsNoImpr = 0;
std::cout << "we enter the initializePaths function" << std::endl;
    initializePaths(allNodes);
std::cout << "we exit the initializePaths function" << std::endl;
check();
std::cout << "we enter the calculateRecord function" << std::endl;
    double record = calculateRecord(allPaths);
std::cout << "we exit the calculateRecord function" << std::endl;
check();
    double prevRecord = record;
    deviation = (1.0-(prob/100))*record;
    for(int k = 0; k < kLim; k++) {
        for(int i = 0; i < iLim; i++) {
            int actiondone = 0;
std::cout << "we enter the twoPointExchange function" << std::endl;
            twoPointExchange();
std::cout << "we exit the twoPointExchange function" << std::endl;
check();
std::cout << "we enter the onePointMove function" << std::endl;
            onePointMove();
std::cout << "we exit the onePointMove function" << std::endl;
check();
std::cout << "we enter the cleanUp function" << std::endl;
            cleanUp();
std::cout << "we exit the cleanUp function" << std::endl;
check();
            if (prevRecord == record) {
                break;
            }
        }
        if (record <= prevRecord) {
            stepsNoImpr++;
        }
        else {
            stepsNoImpr = 0;
        }
        if (stepsNoImpr >= 5) {
            break;
        }
std::cout << "we enter the reinitializationI function" << std::endl;
        reinitializationI(k);
std::cout << "we exit the reinitializationI function" << std::endl;
check();
    }
    p = 2.5;
std::cout << "we enter the reinitializationII function" << std::endl;
    reinitializationII(kLim);
std::cout << "we exit the reinitializationII function" << std::endl;
check();
    for(int k=0; k<kLim; k++) {
        for(int i=0; i<iLim; i++) {
            int actiondone = 0;
            twoPointExchange();
            onePointMove();
            cleanUp();
            if (prevRecord == record) break;
        }
        if (record <= prevRecord) {
            stepsNoImpr++;
        }
        else {
            stepsNoImpr = 0;
        }
        if (stepsNoImpr >= 5) break;
        reinitializationI(k);
    }
    //return best;
check();
}