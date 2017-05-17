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
int bbb;
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
        int id; //.at(0-M) value equal to the cell in the cost array that refers to this node
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

//this function returns the energy cost to visit the node
double visitCost(Node* a) {
    return (a->model.stdDeviation*2)/(sampleRadius*2);
}

//this function takes a pointer to a node and returns the cost to reach it and visit it from the start point
double travelCost(Node* a) {
    double x = a->x;
    double y = a->y;
    return sqrt(x * x + y * y);
}

//this function takes a pointer to a start node "a" and a pointer to a destination node "b" returns the cost to reach it and visit it from the start node
double travelCost(Node* a, Node* b) {
    double x = a->x - b->x;
    double y = a->y - b->y;
    return sqrt(x * x + y * y);
}

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

void check() {
    std::cout << "Result:" << std::endl;
    int flag;
    //std::cout << "allPaths.size(): " << allPaths.size() << std::endl;
    for (int i = 0; i < allPaths.size(); i++) {//FOR EACH PATH
        if (allPaths.at(i)->energySpent > energy) {
            std::cout << "ERRORE" << std::endl;
            flag = 1;
        }
        std::cout << "path " << allPaths.at(i)->id << ": with energySpent: " << allPaths.at(i)->energySpent << "; ";
        for (std::list<Node*>::iterator it = allPaths.at(i)->nodes.begin(); it != allPaths.at(i)->nodes.end(); it++) {//FOR EACH NODE IN THAT PATH
            std::cout << (*it)->id << ": ";
            for (int j = 0; j < (*it)->parts.pathAssignments.size(); j++) {
                if ((*it)->parts.pathAssignments.at(j).first == allPaths.at(i)) {
                    std::cout << "(" << (*it)->parts.pathAssignments.at(j).second.first << "-" << (*it)->parts.pathAssignments.at(j).second.second << ")";
                    std::cout << " ";
                    break;
                }
                else {
                    std::cout << "rilevato problema: " << (*it)->parts.pathAssignments.at(j).first->id << "-" << allPaths.at(i)->id;
                }
                std::cout << " ";
            }
        }
        std::cout << std::endl;
    }
    if (flag == 1) std::cin >> aaa;
}

//TODO: interrogarsi se non sia possibile avere più posizioni ottime di cui una sia preferibile a causa della possibilità di replanning.
//this function takes a pointer to a node, a pointer to a path and returns a pair composed by the minimal cost of inserting the cost into the path and the optimal insertion position.
std::pair<double, int> fitMove(Node* a, Path* path, double visitCost) {
    double minCost = DBL_MAX;
    double cost = DBL_MAX;
    int position = 0;
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
        cost = travelCost(a) + travelCosts.at(a->id).at((*it)->id) - travelCost(*it) + visitCost;
        if (cost + path->energySpent < energy && cost + path->energySpent < energy) {
            minCost = cost;
            position = actual;
            admissible = 1;
        }
        while (itNext != path->nodes.end()) {//POSITION 1 TO LAST
            actual++;
            cost = travelCosts.at((*it)->id).at(a->id) + travelCosts.at(a->id).at((*itNext)->id) - travelCosts.at((*it)->id).at((*itNext)->id) + visitCost;
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
        actual++;
        cost = travelCosts.at((*it)->id).at(a->id) + travelCost(a) - travelCost(*it) + visitCost;
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
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths, std::pair<Path*, std::pair<double, double>>* part) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = 0;
    int next = -1;
    int samePath = 0;;
    double cost;
    int admissible;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//CHECK ON EACH POSITION
            admissible = 0;
            samePath = 0;
            for (std::vector<std::pair<Path*, std::pair<double, double>>>::iterator it = nodes.at(i)->parts.pathAssignments.begin(); it != nodes.at(i)->parts.pathAssignments.end(); it++) {
                if ((*it).first == paths.at(j)) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {
                std::pair<double, int> result = fitMove(nodes.at(i), paths.at(j), part->second.first);//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if (cost < minCost && paths.at(j)->energySpent + cost < energy) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths.at(j)).id;
                        minNode = i;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths.at(j)->id;
                        minNode = nodes.at(i)->id;
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
        std::list<Node*>::iterator it = find(originPath->nodes.begin(), originPath->nodes.end(), allNodes.at(nodeId));
        originPath->totalScore = originPath->totalScore - part->second.second;
            std::list<Node*>::iterator itNext = it;
            itNext++;
        if (it != originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCosts.at((*itPrev)->id).at(nodeId) - travelCosts.at(nodeId).at((*itNext)->id) + travelCosts.at((*itPrev)->id).at((*itNext)->id) - part->second.first;
        }
        if (it == originPath->nodes.begin() && itNext != originPath->nodes.end()) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            originPath->energySpent = originPath->energySpent - travelCost(allNodes.at(nodeId)) - travelCosts.at(nodeId).at((*itNext)->id) + travelCost(*itNext) - part->second.first;
        }
        if (it != originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            originPath->energySpent = originPath->energySpent - travelCosts.at((*itPrev)->id).at(nodeId) - travelCost(allNodes.at(nodeId)) + travelCost(*itPrev) - part->second.first;
        }
        if (it == originPath->nodes.begin() && itNext == originPath->nodes.end()) {
            //THIS LINE OF CODE IS SUBJECT TO ROUNDING ERRORS, BUT WE CAN JUST SET THE ENERGY SPENT TO ZERO
            //originPath->energySpent = originPath->energySpent - (travelCost(allNodes.at(nodeId))*2);
            originPath->energySpent = 0.0;
        }
        originPath->nodes.erase(std::remove(originPath->nodes.begin(), originPath->nodes.end(), allNodes.at(nodeId)), originPath->nodes.end());
    }
    //ADD PART TO THE NEW PATH
    if (destPath >= 0) {//IF THE id OF THE DESTINATION PATH IS POSITIVE THE DESTINATION PATH IS ASSUMED TO EXIST
        std::list<Node*>::iterator it = allPaths.at(destPath)->nodes.begin();
        std::advance(it, position);
        part->first = allPaths.at(destPath);
        if (std::distance(allPaths.at(destPath)->nodes.begin(), it) == allPaths.at(destPath)->nodes.size()) {
            allPaths.at(destPath)->nodes.push_back(allNodes.at(nodeId));
        }
        else {
            allPaths.at(destPath)->nodes.insert(it, allNodes.at(nodeId));
        }
        it = allPaths.at(destPath)->nodes.begin();
        std::advance(it, position);
        allPaths.at(destPath)->totalScore = allPaths.at(destPath)->totalScore + part->second.second;
        if (position != 0 && position != allPaths.at(destPath)->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths.at(destPath)->energySpent = allPaths.at(destPath)->energySpent + travelCosts.at((*itPrev)->id).at(nodeId) + travelCosts.at(nodeId).at((*itNext)->id) - travelCosts.at((*itPrev)->id).at((*itNext)->id) + part->second.first;
        }
        if (position == 0 && position != allPaths.at(destPath)->nodes.size()-1) {
            std::list<Node*>::iterator itNext = it;
            itNext++;
            allPaths.at(destPath)->energySpent = allPaths.at(destPath)->energySpent + travelCost(allNodes.at(nodeId)) + travelCosts.at(nodeId).at((*itNext)->id) - travelCost(*itNext) + part->second.first;
        }
        if (position != 0 && position == allPaths.at(destPath)->nodes.size()-1) {
            std::list<Node*>::iterator itPrev = it;
            itPrev--;
            allPaths.at(destPath)->energySpent = allPaths.at(destPath)->energySpent + travelCosts.at((*itPrev)->id).at(nodeId) + travelCost(allNodes.at(nodeId)) - travelCost(*itPrev) + part->second.first;
        }
        if (position == 0 && position == allPaths.at(destPath)->nodes.size()-1) {
            allPaths.at(destPath)->energySpent = allPaths.at(destPath)->energySpent + travelCost(allNodes.at(nodeId))*2 + part->second.first;
        }
    }
}

//this function takes a vector of pointers to nodes and returns a vector matrix representing the travel+visit costs ordered by node ID number
std::vector<std::vector<double> > updateCosts(std::vector<Node *> allNodes) {
    std::vector<Node *>::iterator it = allNodes.begin();
    //resize the matrix accordingly to the number of nodes
    std::vector<std::vector<double> > distance (allNodes.size());
    for (int i = 0; i < allNodes.size(); i++) {
        distance.at(i).resize(allNodes.size());
    }
    //update the matrix
    for (int i = 0; i < allNodes.size(); i++) {
        for (int j = 0; j < allNodes.size(); j++) {
            distance.at(i).at(j) = travelCost(allNodes.at(i), allNodes.at(j));
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
            //the subsequent lines are of the form .at(x, y, score) where x is the position on the x-axis, y is the position on the y-axis and score is the score associated with the target
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
            if ((travelCost(newNode) + visitCost(newNode)) < energy) {
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

//this function takes a vector of nodes, a vector of paths and returns the cheapest insertion of a node from the first vector to a path in the second in the form of pair of pairs. The first pair is the chosen node's id and the cost. The second one is the path and the position where to insert.
std::pair< std::pair<int, double>, std::pair<int, int> > cheapestInsertion(std::vector<Node*> nodes, std::vector<Path*> paths) {
    double minCost = DBL_MAX;
    int minPath = -1;
    int minNode = -1;
    int position = 0;
    int next = -1;
    double cost;
    int admissible;
    int samePath;
    for (int i = 0; i < nodes.size(); i++) {//FOR EACH NODE
        for (int j = 0; j < paths.size(); j++) {//FOR EACH PATH
            admissible = 0;
            samePath = 0;
            for (int l = 0; l < nodes.at(i)->parts.pathAssignments.size(); l++) {//FOR EACH PART BELONGING TO THE NODE
                if ((nodes.at(i)->parts.pathAssignments.at(l).first == paths.at(j)) && (nodes.at(i)->parts.pathAssignments.at(l).first != NULL)) {
                    samePath = 1;
                    break;
                }
            }
            if (samePath != 1) {//IF THE PART IS NOT ALREADY IN THE NODE
                std::pair<double, int> result = fitMove(nodes.at(i), paths.at(j), visitCost(nodes.at(i)));//USE fitMove TO GET FITNESS OF CHEAPEST INSERTION OF NODE CONSIDERED IN PATH CONSIDERED
                cost = result.first;
                if (admissible == 0) {//IF SOLUTION IS ADMISSIBLE AND BETTER THAN PREVIOUS WE SAVE IT
                    if ((cost < minCost) && ((paths.at(j)->energySpent + cost) < energy)) {//WE PERFORM THE CHECK OF ADMISSIBILITY ONLY IF WE DIDN'T ALREADY FIND AN ADMISSIBLE SOLUTION FOR THIS PATH
                        position = result.second;
                        minCost = cost;
                        minPath = (*paths.at(j)).id;
                        minNode = nodes.at(i)->id;
                        admissible = 1;
                    }
                }
                else {
                    if (cost < minCost) {
                        position = result.second;
                        minCost = cost;
                        minPath = paths.at(j)->id;
                        minNode = nodes.at(i)->id;
                    }
                }
            }
        }
    }
    return std::make_pair(std::make_pair(minNode, minCost), std::make_pair(minPath, position));
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
        sum = sum + tmpAllPaths.at(i)->totalScore;
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
        nodesToAdd.at(i)->parts.pathAssignments.push_back(std::make_pair(pathPtr, std::make_pair(visitCost(nodesToAdd.at(i)), nodesToAdd.at(i)->score)));
        nodesToAdd.at(i)->parts.pathAssignments.at(0).first = NULL;
        movePart(&nodesToAdd.at(i)->parts.pathAssignments.at(0), nodesToAdd.at(i)->id, pathPtr->id, 0);
        nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), nodesToAdd.at(i)), nodesToAdd.end());
    }
    int prevSize;
    while (nodesToAdd.size() > 0) {//WHILE THERE ARE NODES WITHOUT PATH WE INSERT ONE OF THOSE IN A PATH, WHEN NECESSARY WE CREATE A NEW ONE
        prevSize = nodesToAdd.size();
        std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(nodesToAdd, allPaths);
        if (result.first.first == -1) {//IF THERE'S NO ADMISSIBLE INSERTION THEN CREATE A NEW PATH AND REPEAT
            pathPtr = new Path(allPaths.size());
        }
        else {//ELSE ADD THE NODE TO THE PATH, UPDATE THE NODE'S PATH FIELD AND THE PATH'S SCORE AND ENERGY SPENT FIELDS
            std::pair< Path*, std::pair<double, double> > tmp = std::make_pair(allPaths.at(result.second.first), std::make_pair(visitCost(allNodes.at(result.first.first)), allNodes.at(result.first.first)->score));
            allNodes.at(result.first.first)->parts.pathAssignments.push_back(tmp);
            allNodes.at(result.first.first)->parts.pathAssignments.at(0).first = NULL;
            movePart(&allNodes.at(result.first.first)->parts.pathAssignments.at(0), result.first.first, result.second.first, result.second.second);
            nodesToAdd.erase(std::remove(nodesToAdd.begin(), nodesToAdd.end(), allNodes.at(result.first.first)), nodesToAdd.end());
        }
    }
    record = calculateRecord(allPaths);//SOME NON-TRIVIAL WORK WAS DONE ON THE allPaths VECTOR SO WE UPDATE THE record
}

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
        for (int l = 0; (l < allNodes.at(i)->parts.pathAssignments.size()) && (actionDone == 0); l++) {//FOR EACH PART
            candPart = &allNodes.at(i)->parts.pathAssignments.at(l);
            prevPath = candPart->first->id;
            prevPosition = 0;
            for (std::list<Node*>::iterator tmpIt = allPaths.at(prevPath)->nodes.begin(); tmpIt != allPaths.at(prevPath)->nodes.end(); tmpIt++) {
                if (*(tmpIt) == allNodes.at(i)) {//IF WE'VE FOUND THE POSITION WE EXIT THE FOR LOOP
                    break;
                }
                prevPosition++;
            }
            for (int j = 0; (j < allPaths.size()) && (actionDone == 0); j++) {//WE TRY TO MOVE IT IN ANOTHER PATH
                Path* candPath = allPaths.at(j);
                if (candPath->id != prevPath) {
                    std::pair<double, int> result = fitMove(allNodes.at(i), candPath, candPart->second.first);
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
                                bestMovePart = &allNodes.at(i)->parts.pathAssignments.at(l);
                            }
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
bool comparePartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>*> a, std::pair<int, std::pair<Path*, std::pair<double, double>>*> b) {
    if (a.second->second.second < b.second->second.second) {
        return false;
    }
    else {
        return true;
    }
}

//COMPARISON FUNCTION BETWEEN PARTS BASED ON SCORE WEIGHTED ON THE LOCAL COST OF THAT PART IN THE PATH (INCLUDING VISITCOST)
bool compareWeightedPartScore(std::pair<int, std::pair<Path*, std::pair<double, double>>*> a, std::pair<int, std::pair<Path*, std::pair<double, double>>*> b) {
    std::list<Node*>::iterator itA =  find(a.second->first->nodes.begin(), a.second->first->nodes.end(), allNodes.at(a.first));
    std::list<Node*>::iterator itALast = a.second->first->nodes.end();
    itALast--;
    std::list<Node*>::iterator itB =  find(b.second->first->nodes.begin(), b.second->first->nodes.end(), allNodes.at(b.first));
    std::list<Node*>::iterator itBLast = b.second->first->nodes.end();
    itBLast--;
    double costA = 0.0;
    if ((itA != a.second->first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        aNext++;
        double costA = travelCosts.at((*aPrev)->id).at((*itA)->id) + travelCosts.at((*itA)->id).at((*aNext)->id) - travelCosts.at((*aPrev)->id).at((*aNext)->id) + a.second->second.first;
    }
    if ((itA != a.second->first->nodes.begin()) && (itA == itALast)) {
        std::list<Node*>::iterator aPrev = itA;
        aPrev--;
        double costA = travelCosts.at((*aPrev)->id).at((*itA)->id) + travelCost(allNodes.at((*itA)->id)) - travelCost(allNodes.at((*aPrev)->id)) + a.second->second.first;
    }
    if ((itA == a.second->first->nodes.begin()) && (itA != itALast)) {
        std::list<Node*>::iterator aNext = itA;
        aNext++;
        double costA = travelCost(allNodes.at((*itA)->id)) + travelCosts.at((*itA)->id).at((*aNext)->id) - travelCost(allNodes.at((*aNext)->id)) + a.second->second.first;
    }
    if ((itA == a.second->first->nodes.begin()) && (itA == itALast)) {
        double costA = travelCost(allNodes.at((*itA)->id)) + travelCost(allNodes.at((*itA)->id)) + a.second->second.first;
    }
    double costB = 0.0;
    if ((itB != b.second->first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        bNext++;
        double costB = travelCosts.at((*bPrev)->id).at((*itB)->id) + travelCosts.at((*itB)->id).at((*bNext)->id) - travelCosts.at((*bPrev)->id).at((*bNext)->id) + b.second->second.first;
    }
    if ((itB != b.second->first->nodes.begin()) && (itB == itBLast)) {
        std::list<Node*>::iterator bPrev = itB;
        bPrev--;
        double costB = travelCosts.at((*bPrev)->id).at((*itB)->id) + travelCost(allNodes.at((*itB)->id)) - travelCost(allNodes.at((*bPrev)->id)) + b.second->second.first;
    }
    if ((itB == b.second->first->nodes.begin()) && (itB != itBLast)) {
        std::list<Node*>::iterator bNext = itB;
        bNext++;
        double costB = travelCost(allNodes.at((*itB)->id)) + travelCosts.at((*itB)->id).at((*bNext)->id) - travelCost(allNodes.at((*bNext)->id)) + b.second->second.first;
    }
    if ((itB == b.second->first->nodes.begin()) && (itB == itBLast)) {
        double costB = travelCost(allNodes.at((*itB)->id)) + travelCost(allNodes.at((*itB)->id)) + b.second->second.first;
    }
    return (a.second->second.second/costA < b.second->second.second/costB);
}

//THIS FUNCTION IMPLEMENTS THE CLEAN UP ROUTINE
void cleanUp() {
    double record = 0.0;
    int actionDone = 0;
    double cost = 0.0;
    for (int p = 0; p < allPaths.size(); p++) { //2-OPT IS APPLIED FOR EVERY PATH
        int pathNodeNumber = allPaths.at(p)->nodes.size();
        while (cost < 0.0) {//WHILE COST LESS THAN ZERO THERE'S AN IMPROVEMENT IN ENERGY COST
            cost = 0.0;
                std::list<Node*>::iterator firstNode = allPaths.at(p)->nodes.begin();
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
                            cost = cost - travelCosts.at((*firstNodePrev)->id).at((*firstNode)->id) + travelCosts.at((*firstNodePrev)->id).at((*secondNode)->id);
                        }
                        if (j == pathNodeNumber-1) {
                            cost = cost - sqrt((*secondNode)->x*(*secondNode)->x + (*secondNode)->y*(*secondNode)->y) + sqrt((*firstNode)->x*(*firstNode)->x + (*firstNode)->y*(*firstNode)->y);
                        }
                        else {
                            std::list<Node*>::iterator secondNodeNext = secondNode;
                            secondNodeNext++;
                            cost = cost - travelCosts.at((*secondNode)->id).at((*secondNodeNext)->id) + travelCosts.at((*firstNode)->id).at((*secondNodeNext)->id);
                        }
                        std::list<Node*>::iterator secondNodePrev = secondNode;
                        secondNodePrev--;
                        std::list<Node*>::iterator firstNodeNext = firstNode;
                        firstNodeNext++;
                        cost = cost - travelCosts.at((*firstNode)->id).at((*firstNodeNext)->id) + travelCosts.at((*secondNode)->id).at((*firstNodeNext)->id);
                        cost = cost - travelCosts.at((*secondNodePrev)->id).at((*secondNode)->id) + travelCosts.at((*secondNodePrev)->id).at((*firstNode)->id);
                        if (cost < 0.0) {
                            std::list<Node*> copy;
                            copy.assign(firstNode, secondNode);
                            copy.reverse();
                            allPaths.at(i)->nodes.erase(firstNode, secondNode);
                            allPaths.at(i)->nodes.insert(firstNode, copy.begin(), copy.end());
                            allPaths.at(i)->energySpent = allPaths.at(i)->energySpent + cost;
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
        for (std::list<Node*>::iterator it = tmpAllPaths.at(i)->nodes.begin(); it != tmpAllPaths.at(i)->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes.at(i)->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes.at(i)->id, &topNodes.at(i)->parts.pathAssignments.at(j)));
        }
    }
    std::sort(topParts.begin(), topParts.end(), comparePartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts.at(i).first < allNodes.size()); i++) {
        dropPart(allNodes.at(topParts.at(i).first), topParts.at(i).second);
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
        for (std::list<Node*>::iterator it = tmpAllPaths.at(i)->nodes.begin(); it != tmpAllPaths.at(i)->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes.at(i)->parts.pathAssignments.size(); j++) {
            topParts.push_back(std::make_pair(topNodes.at(i)->id, &topNodes.at(i)->parts.pathAssignments.at(j)));
        }
    }
    std::sort(topParts.begin(), topParts.end(), compareWeightedPartScore);
    for (int i = 0; (i < k) && (i < topParts.size()) && (topParts.at(i).first < allNodes.size()); i++) {
        dropPart(allNodes.at(topParts.at(i).first), topParts.at(i).second);
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
    int candPart;
    //CREATE A VECTOR OF THE NODES IN TOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> topNodes;
    //WE CREATE A SORTED VERSION OF allPaths
    std::vector<Path*> tmpAllPaths = allPaths;
    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
    //CREATE A VECTOR OF REFERENCES TO THE NODES IN THE FIRST M PATH (THE TOP ONES)
    for (int i = 0; i < M; i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths.at(i)->nodes.begin(); it != tmpAllPaths.at(i)->nodes.end(); it++) {
            topNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS OF THE TOP NODES, SAVED IN THE topNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> topParts;
    for (int i = 0; i < topNodes.size(); i++) {
        for (int j = 0; j < topNodes.at(i)->parts.pathAssignments.size(); j++) {
            std::pair <int, std::pair<Path*, std::pair<double, double>>*> tmp = std::make_pair(topNodes.at(i)->id, &(topNodes.at(i)->parts.pathAssignments.at(j)));
            topParts.push_back(tmp);
        }
    }
    //CREATE A VECTOR OF THE NODES IN NTOP PATHS, SAVED IN THE PATH ORDER
    std::vector<Node*> nTopNodes;
    for (int i = M; i < tmpAllPaths.size(); i++) {
        for (std::list<Node*>::iterator it = tmpAllPaths.at(i)->nodes.begin(); it != tmpAllPaths.at(i)->nodes.end(); it++) {
            nTopNodes.push_back(*it);
        }
    }
    //CREATE A VECTOR OF REFERENCES TO THE PARTS IN NTOP PATHS, SAVED IN THE nTopNodes ORDER
    std::vector<std::pair<int, std::pair<Path*, std::pair<double, double>>*>> nTopParts;
    for (int i = 0; i < nTopNodes.size(); i++) {
        for (int j = 0; j < nTopNodes.at(i)->parts.pathAssignments.size(); j++) {
            std::pair <int, std::pair<Path*, std::pair<double, double>>*> tmp = std::make_pair(nTopNodes.at(i)->id, &(nTopNodes.at(i)->parts.pathAssignments.at(j)));
            nTopParts.push_back(tmp);
        }
    }
    for (int i = 0; i < topParts.size(); i++) {//FOR EACH PART IN topParts
        candScore = -DBL_MAX;
        candPart = -1;
        candNode = -1;
        candDestPath = -1;
        candPositionForTopNode = -1;
        candPositionForNTopNode = -1;
        actionDone = 0;
        totScore = -DBL_MAX;
        //REMOVE THE PART FROM THE TOP PATH IN ORDER TO GET AN ACCURATE ESTIMATE OF THE COST BY USING fitMove and cheapestInsertion
        Path* topOrigPath = topParts.at(i).second->first;
        std::list<Node*>::iterator topOrigPositionIt = find(topOrigPath->nodes.begin(), topOrigPath->nodes.end(), allNodes.at(topParts.at(i).first));
        int topOrigPosition = distance(topOrigPath->nodes.begin(), topOrigPositionIt);
        movePart(topParts.at(i).second, topParts.at(i).first, -1, 0);
        for (int j = 0; (j < nTopParts.size()) && (actionDone == 0); j++) {//FOR EACH PART IN nTopParts
            if (topParts.at(i).first != nTopParts.at(j).first) {//CHECK THAT THE TWO PARTS ARE NOT OF THE SAME NODE
                //CHECK IF THE EXCHANGE IS POSSIBLE, NOTE THAT WE NEED TO CHECK ONLY THE TOP NODE SIDE BECAUSE WE CAN ALWAYS CREATE A NEW NTOP PATH
                std::pair<double, int> toTop = fitMove(allNodes.at(nTopParts.at(j).first), topOrigPath, nTopParts.at(j).second->second.first);
                cost = toTop.first;
                topPathPos = toTop.second;
                if (cost + topOrigPath->energySpent < energy) {//IF THE NTOP PART FITS WE PROCEED TO PERFORM THE TWO movePart
                    //MOVE NTOP PART IN TOP PATH
                    Path* nTopOrigPath = nTopParts.at(j).second->first;
                    std::list<Node*>::iterator nTopOrigPositionIt = find(nTopOrigPath->nodes.begin(), nTopOrigPath->nodes.end(), allNodes.at(nTopParts.at(j).first));
                    int nTopOrigPosition = distance(nTopOrigPath->nodes.begin(), nTopOrigPositionIt);
                    movePart(nTopParts.at(j).second, nTopParts.at(j).first, topOrigPath->id, topPathPos);
                    //MOVE TOP PART IN NTOP PATH
                    std::vector<Path*> tmpAllPaths = allPaths;
                    std::sort(tmpAllPaths.begin(), tmpAllPaths.end(), scorePathComp);
                    std::vector<Path*>::iterator firstNTopPath = tmpAllPaths.begin();
                    std::advance(firstNTopPath, M);
                    std::vector<Path*> nTopPaths (firstNTopPath, tmpAllPaths.end());
                    std::pair< std::pair<int, double>, std::pair<int, int> > result = cheapestInsertion(allNodes.at(topParts.at(i).first), nTopPaths, topParts.at(i).second);
                    if (result.second.first != -1) {
                        movePart(topParts.at(i).second, topParts.at(i).first, result.second.first, result.second.second);
                    }
                    else {
                        Path* newNTopPath = new Path(allPaths.size());
                        movePart(topParts.at(i).second, topParts.at(i).first, newNTopPath->id, result.second.second);
                    }
                    totScore = calculateRecord(allPaths);
                    if (totScore > record) {
                        actionDone = 1;
                        nTopParts.erase(std::remove(nTopParts.begin(), nTopParts.end(), nTopParts.at(j)), nTopParts.end());
                    }
                    else {
                        if (totScore > candScore) {//SAVE BEST MOVE
                            candScore = totScore;
                            candPart = j;
                            candDestPath = result.second.first;
                            candPositionForTopNode = result.second.second;
                            candPositionForNTopNode = topPathPos;
                        }
                        //REVERT THE TWO movePart
                        movePart(topParts.at(i).second, topParts.at(i).first, -1, 0);
                        movePart(nTopParts.at(j).second, nTopParts.at(j).first, nTopOrigPath->id, nTopOrigPosition);
                    }
                }
            }
        }
        if (actionDone == 0) {
            if (candScore > deviation) {
                movePart(nTopParts.at(candPart).second, nTopParts.at(candPart).first, topOrigPath->id, candPositionForNTopNode);
                if (candDestPath == -1) {
                    Path* tmpPath = new Path(allPaths.size());
                    movePart(topParts.at(i).second, topParts.at(i).first, tmpPath->id, 0);
                }
                else {
                    movePart(topParts.at(i).second, topParts.at(i).first, candDestPath, candPositionForTopNode);
                }
                nTopParts.erase(std::remove(nTopParts.begin(), nTopParts.end(), nTopParts.at(candPart)), nTopParts.end());
            }
            else {
                movePart(topParts.at(i).second, topParts.at(i).first, topOrigPath->id, topOrigPosition);
            }
        }
    }
}

int main(int argc, char** argv) {
    initializeNodes(argv[1]);
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
            if (prevRecord == record) break;
        }
        if (record <= prevRecord) {
            stepsNoImpr++;
        }
        else {
            stepsNoImpr = 0;
        }
        if (stepsNoImpr >= 5) break;
std::cout << "we enter the reinitializationI function" << std::endl;
        reinitializationI(k);
std::cout << "we exit the reinitializationI function" << std::endl;
check();
    }
    //TODO write best solution to file;
std::cout << "FINAL" << std::endl;
check();
}
