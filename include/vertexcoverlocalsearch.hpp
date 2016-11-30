//
//  vertexcoverlocalsearch.hpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 01.06.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#ifndef vertexcoverlocalsearch_hpp
#define vertexcoverlocalsearch_hpp

#include <stdio.h>
#include "mkgraph.hpp"

using namespace std;

/**
 *  Default value that is used to count how often the second pruning rule is
 *  applied. The corresponding parameter is passed by reference, which is why
 *  its default value has to be declared statically here.
 */
static long long pruningRule1NumberOfApplicationsDefault = 0;
static long long pruningRule2NumberOfApplicationsDefault = 0;
static long long pruningRule3NumberOfApplicationsDefault = 0;
static long long isolatedVerticesFoundWhileApplyingPruningRule3Default = 0;

class VertexCoverLocalSearch
{
private:
    
#pragma mark - Tree Structure
    
    /**
     *  We want to keep track of how many black / white descendants a vertex in our
     *  exchange set has. This will enable us to use some pruning rules that
     *  consider properties of certain subtrees in the exchange set.
     */
    struct MKVertexDescendant
    {
        Graph::vertex vertex;       // Vertex itself.
        Graph::vertex parent;       // Parent of the vertex in the tree.
        int blackDescendants = 0;   // Number of black descendants of the vertex in its subtree.
        int whiteDescendants = 0;   // Number of white descendants of the vertex in its subtree.
        
        void addBlackDescendant(map<Graph::vertex, MKVertexDescendant> &descendants)
        {
            blackDescendants += 1;
            
            /**
             *  Now update all its ancestors values.
             */
            if (vertex != parent)
            {
                descendants[parent].addBlackDescendant(descendants);
            }
        }
        
        void addWhiteDescendant(map<Graph::vertex, MKVertexDescendant> &descendants)
        {
            whiteDescendants += 1;
            
            /**
             *  Now update all its ancestors values.
             */
            if (vertex != parent)
            {
                descendants[parent].addWhiteDescendant(descendants);
            }
        }
    };
    
    static void addEdgesBetweenVerticesToGraph(vector<pair<Graph::vertex, Graph::vertex_set>> &vertices, MKGraph &G);
    
    /**
     *  Method that is used to add a black vertex to the tree greedily. That means
     *  v is a vertex that can be added without having to add new white vertices to
     *  the tree.
     *
     *  @param v                 Vertex to add greedily
     *  @param G                 Graph that we're searching a vertex cover in
     *  @param C                 Known vertex cover
     *  @param T                 Current Tree
     *  @param B                 Set of black vertices that must not be processed further.
     *  @param TB                Set of black vertices in T.
     *  @param P                 Set of potential black vertices that could still be added to T.
     *  @param descendants       Data structure that counts the number of black / white descendants that a vertex has.
     *  @param applyPruningRule2 Flag that describes whether the second pruning rule should be applied.
     */
    static void addBlackVertexToTreeGreedily(Graph::vertex &v,
                                             const MKGraph &G,
                                             const Graph::vertex_set &C,
                                             Graph::vertex_set &T,
                                             Graph::vertex_set &B,
                                             Graph::vertex_set &TB,
                                             Graph::vertex_set &P,
                                             map<Graph::vertex, MKVertexDescendant> &descendants,
                                             bool applyPruningRule2);
    
    /**
     *  Sorts the vertices of a subset of the vertices in a graph by their degree, descending.
     *
     *  @param G               Graph that contains the vertices and the vertex set.
     *  @param vertexSet       Vertex set to sort by degree.
     *  @param sortedVertexSet Array that will contain the vertices of the vertex set sorted by degree, descending.
     */
    static void sortVerticesFromSetInGraphByDegree(vector<int> &vertexDegrees,
                                                   const Graph::vertex_set &vertexSet,
                                                   vector<Graph::vertex> &sortedVertexSet);
    
    static Graph::vertex_set blackVerticesFromSet(const Graph::vertex_set &C,
                                                  const std::vector<bool> &vertexColors);
    
    static Graph::vertex_set whiteVerticesFromSet(const Graph::vertex_set &C,
                                                  const std::vector<bool> &vertexColors);
    
    /**
     *  This method is the equivalent to the ENUM method in the paper. It is used
     *  to recursively enumerate potential exchange sets when searching for a better
     *  vertex cover.
     *
     *  @param G                                              Graph in which the vertex sets are to be enumerated.
     *  @param C                                              Known vertex cover that is to be improved.
     *  @param k                                              Number of allowed exchanges, which is also the size of the sets that we enumerate.
     *  @param u                                              Pivot vertex that is to be processed now.
     *  @param T                                              Vertex set that we're currently enumerating.
     *  @param B                                              Set containing black vertices that must not be processed anymore.
     *  @param W                                              Set of white vertices that must not be processed anymore.
     *  @param TB                                             Set of black vertices in T.
     *  @param TW                                             Set of white vertices in T.
     *  @param P                                              Set of potential black vertices that could still be added to T.
     *  @param pivotStack                                     Stack containing the white pivot vertices that still have to be processed.
     *  @param searchNodeCount                                Used to keep track of how often this method is called recursively.
     *  @param descendants                                    To allow for certain pruning rules to be applied we want to know for each vertex how many black and white descendants it has.
     *  @param applyPruningRule1                              Bool that decides whether the first pruning rule should be applied. Default is false.
     *  @param pruningRule1NumberOfApplications               Counts how often the first pruning rule was applied. Default value is 0. Will not be increased if applyPruningRule1 = false.
     *  @param applyPruningRule2                              Bool that decides whether the second pruning rule should be applied. Default is false. (Subtrees with negative contributions)
     *  @param pruningRule2NumberOfApplications               Counts how often the second pruning rule was applied. Default value is 0. Will not be increased if applyPruningRule2 = false.
     *  @param applyPruningRule3                              If set to true the third pruning rule (Greedy Matching) will be applied. Default is false.
     *  @param pruningRule3NumberOfApplications               Counts how often the third pruning rule (Greedy Matching) was applied. Default value is 0. Will not be increased if usePruningRule3 = false.
     *  @param isolatedVerticesFoundWhileApplyingPruningRule3 Counts how often the isolated vertices were found while applying the third pruning rule 0. Will not be increased if usePruningRule3 = false.
     *
     *  @return True, when the found set T is a valid exchange set. False, else.
     */
    static bool enumerateKExchangeSets(const MKGraph &G,
                                       const Graph::vertex_set &C,
                                       const std::vector<bool> &vertexColors,
                                       int k,
                                       Graph::vertex u,
                                       Graph::vertex_set &T,
                                       Graph::vertex_set &B,
                                       Graph::vertex_set &W,
                                       Graph::vertex_set &TB,
                                       Graph::vertex_set &TW,
                                       Graph::vertex_set &P,
                                       stack<Graph::vertex> &pivotStack,
                                       long long &searchNodeCount,
                                       map<Graph::vertex, MKVertexDescendant> &descendants,
                                       bool applyPruningRule1,
                                       long long &pruningRule1NumberOfApplications,
                                       bool applyPruningRule2,
                                       long long &pruningRule2NumberOfApplications,
                                       bool applyPruningRule3,
                                       long long &pruningRule3NumberOfApplications,
                                       long long &isolatedVerticesFoundWhileApplyingPruningRule3,
                                       long long &trianglesFoundWhileApplyingPruningRule3);
    
    /**
     *  Given a vertex Cover C, this method tries to find a vertex cover Cprime
     *  that has one vertex less than C by exchanging at most k vertices.
     *
     *  @param G                            Graph in which the vertex cover is to be improved.
     *  @param C                            Vertex cover to be improved.
     *  @param kMax                         Maximum number of exchanges to be done.
     *  @param Cprime                       Improved solution, or an empty set if a better solution could not be found.
     *  @param elementToStartSearchAt       Defines at which element in the list of cover vertices that are sorted by degree. If -1, we start with the first element in the cover.
     *  @param elementToContinueSearchAt    Will be the element that would have been the next candidate when an improvement was found.
     *  @param applyPruningRule1            If set to true the first pruning rule will be applied. Default is false.
     *  @param applyPruningRule2            If set to true the second pruning rule (Positive contribution of subtree) will be applied. Default is false.
     *  @param applyPruningRule3            If set to true the third pruning rule (White vertex max reached) will be applied. Default is false.
     *
     *  @return true, when a better vertex set Cprime was found. False, else.
     */
    static bool improveVertexCoverLocallyByOne(MKGraph &G,
                                               vector<int> &vertexDegrees,
                                               Graph::vertex_set &C,
                                               int kMax,
                                               Graph::vertex_set &Cprime,
                                               int elementToStartSearchAt,
                                               int &elementToContinueSearchAt,
                                               bool applyPruningRule1 = false,
                                               bool applyPruningRule2 = false,
                                               bool applyPruningRule3 = false,
                                               bool performSanityCheck = false);
    
public:
    
    /**
     *  Calculates a vertex cover by iteratively removing the vertex with
     *  maximum degree from the graph.
     *
     *  @param cover A vertex cover of the graph.
     */
    static void getGreedyVertexCoverFromGraph(const MKGraph &G, MKGraph::vertex_set &cover);
    
    /**
     *  Given a vertex Cover C, this method tries to find a vertex cover Cprime
     *  that is as small as possible, by iteratively trying to reduce the vertex
     *  cover size by 1.
     *
     *  @param G                    Graph in which the vertex cover is to be improved.
     *  @param C                    Vertex cover to be improved.
     *  @param kMax                 Maximum number of exchanges to be done.
     *  @param Cprime               Improved solution, or an empty set if a better solution could not be found.
     *  @param applyPruningRule1    If set to true the first pruning rule will be applied. Default is false.
     *  @param applyPruningRule2    If set to true the second pruning rule (Positive contribution of subtree) will be applied. Default is false.
     *  @param applyPruningRule3    If set to true the third pruning rule (White vertex max reached) will be applied. Default is false.
     *
     *  @return true, when a better vertex set Cprime was found. False, else.
     */
    static bool improveVertexCoverLocally(const MKGraph &G,
                                          Graph::vertex_set &C,
                                          int kMax,
                                          Graph::vertex_set &Cprime,
                                          bool applyPruningRule1 = false,
                                          bool applyPruningRule2 = false,
                                          bool applyPruningRule3 = false,
                                          bool performSanityCheck = false);
};


#endif /* vertexcoverlocalsearch_hpp */