//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <functional>
#include <limits> //For double Inf
#include <numeric> //For std::accumulate
#include <queue>
#include <cmath>
#include "planning/AStarPlanner.h"
using namespace std;

#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 0.5f
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

namespace SteerLib
{
	AStarPlanner::AStarPlanner(){}

	AStarPlanner::~AStarPlanner(){}

	bool AStarPlanner::canBeTraversed ( int id ) 
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x,z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = MAX(x-OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z-OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());


		for (int i = x_range_min; i<=x_range_max; i+=GRID_STEP)
		{
			for (int j = z_range_min; j<=z_range_max; j+=GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords( i, j );
				traversal_cost += gSpatialDatabase->getTraversalCost ( index );
				
			}
		}

		if ( traversal_cost > COLLISION_COST ) 
			return false;
		return true;
	}

	inline float euclideanDist(const Util::Point& p1, const Util::Point& p2)
	{
		return Util::distanceBetween(p1, p2);
	}

	inline float manhattanDist(const Util::Point& p1, const Util::Point& p2)
	{
		return (abs(p2.x - p1.x) + abs(p2.y - p1.y) + abs(p2.z - p1.z));
	}


	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}


	std::vector<Util::Point> AStarPlanner::reconstructPath(AStarPlannerNode* currentNode)
	{
		std::vector<Util::Point> path;

		while(currentNode != NULL)
		{
			path.insert(path.begin(), currentNode->point);
			currentNode = currentNode->parent;
		}

		return path;
	}

	std::vector<Util::Point> AStarPlanner::getNeighbors(const Util::Point& p)
    {
        int start_x = MAX((p.x-1),gSpatialDatabase->getOriginX());
        int end_x = MIN((p.x+1), gSpatialDatabase->getNumCellsX() + gSpatialDatabase->getOriginX());

        int start_z = MAX((p.z-1), gSpatialDatabase->getOriginZ());
        int end_z = MIN((p.z+1), gSpatialDatabase->getNumCellsZ() + gSpatialDatabase->getOriginZ());

        std::vector<Util::Point> neighbors;

        for (int i = start_x; i <= end_x; ++i) {
            for (int j = start_z; j <= end_z; ++j) {
                if (i != p.x || j != p.z) {
                    if (canBeTraversed(gSpatialDatabase->getCellIndexFromLocation(i, j)))
                    {
                  	   neighbors.push_back(Util::Point(i, 0, j));
                    }
                }
            }
        }
        return neighbors;
    }

   AStarPlannerNode* getNodeByLowestFScore(std::vector<AStarPlannerNode*> openset)
	{
		double minF = std::numeric_limits<double>::max();
		//double minG = std::numeric_limits<double>::max();
		double maxG = -1 * std::numeric_limits<double>::max();

		AStarPlannerNode* returnNode;

		for(std::vector<AStarPlannerNode*>::iterator iter = openset.begin(); iter != openset.end(); ++iter)
		{
			AStarPlannerNode* tmp = (*iter);

			if(tmp->f < minF){
				returnNode = tmp;
				minF = tmp->f;
			}
			else if (tmp->f == minF)
			{
				if(tmp->g > maxG) //(tmp->g < minG) 
				{
					returnNode = tmp;
					maxG = tmp->g; //minG = tmp->g;
				}
			}
		}

		return returnNode;
	}


	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path,  Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		gSpatialDatabase = _gSpatialDatabase;

  		int numberOfOpenNodes = 0;
		/*
		 * heuristicFunc is of type pointer to function taking in to const Util::Point refrences and returning a float
		 *
		 * obviously this will be the pointer to either the euclidean or manhattan distance function.
		 * this function pointer will be set here and passed to all functions which need the heuristic
		 */
		float (*heuristic)(const Util::Point&, const Util::Point&) = manhattanDist;
		std::vector<AStarPlannerNode*> openset;
		std::vector<AStarPlannerNode*> closedset;

		/*
		 * The openset is initialized to have the start point in it
		 *
		 * The constructor for the AStarPlannerNode takes four arguments:
		 * The point: In this case it is the start
		 * The 'g value' which is the distance from the start to the point
		 * 		In this case it is 0 as the point in question is the start point
		 * The 'f value' which is the 'g value' plus the heuristic
		 * 		In this case the 'g value' is 0 so this is simply the value 
		 * 		of the heuristic function to the goal
		 * The 'parent node', this is to be able to rewind through the list once the goal is found
		 * 		Since this is the start node the parent is set to null
		 */

		openset.push_back(new AStarPlannerNode(start, 0, heuristic(start, goal), NULL));

		while(!openset.empty()) {
			//Gets the node in the openset with the lowest f score
			//
			//note operator< is overloaded to use f scores
			AStarPlannerNode* curr = getNodeByLowestFScore(openset);

			if(curr->point == goal) {
				agent_path = reconstructPath(curr);

				//Value for output
				int lengthOfSol = (int) agent_path.size();

				return true;
			}

			//Remove the current node from the openset and put it in the closedset
			openset.erase(std::remove(openset.begin(), openset.end(), curr), openset.end());
			closedset.push_back(curr);

			std::vector<Util::Point> neighbors = getNeighbors(curr->point);
			neighbors.push_back(curr->point);

			for(std::vector<Util::Point>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
				//If the current neighbor is already in the closed set, the best path to it has been found
				//and there is no need to continue working with it
				

				Util::Point currNeighbor = *it;

				if(std::find_if(closedset.begin(), closedset.end(), [&currNeighbor] (AStarPlannerNode* a) -> bool { return a->point==currNeighbor; }) != closedset.end() )
					continue;

				/*
				 * if the openset does not contain the point, this is the first time it has been visited
				 * add it to the openset.
				 *
				 * if it has been previously visited, we need to check if the current path to the neighbor is
				 * better than the previous path. The better path has the lower g score.  If the new g score
				 * is lesser than the currently set g score on that node, update it with new values
				 */
				
				float newGScore = curr->g + heuristic(curr->point, currNeighbor);

				std::vector<AStarPlannerNode*>::iterator node = \
					std::find_if(openset.begin(), openset.end(), [&currNeighbor] (AStarPlannerNode* a) -> bool { return a->point==currNeighbor; });

				if(node == openset.end()) {
					numberOfOpenNodes++;
					openset.push_back(new AStarPlannerNode(currNeighbor, newGScore, newGScore + heuristic(currNeighbor, goal), curr));
				} else if((*node)->g > newGScore) {
					(*node)->g = newGScore;
					(*node)->f = newGScore + heuristic(currNeighbor, goal);
					(*node)->parent = curr;
				}
				
			}

		}

		return false;
	}
}
