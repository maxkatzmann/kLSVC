//
//  vertexcoverlocalsearch.cpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 01.06.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#include "vertexcoverlocalsearch.hpp"

#include <unistd.h>

#pragma mark - Sorting

/**
 *  This struct is used to sort vertices by their degree.
 *
 *  We want to use the STL sort function but need to compare vertices by
 *  degree instead of their ID.
 *
 *  Therefore the MKVertex stores the vertexID and the number of
 *  its neighbors. Its comparator compares MKVertex elements by the number
 *  of their neighbors.
 */
struct MKVertex
{
    /**
     *  Original ID of the vertex.
     */
    unsigned int vertexID;

    /**
     *  The number of its neighbors.
     */
    long neighbors;

    /**
     *  Comparator that compares to MKVertex elements by their neighbors.
     *  Will sort in descending order!
     *
     *  @param lhs MKVertex representing the lefthandside of the comparison operation.
     *  @param rhs MKVertex representing the righthandside of the comparison operation.
     *
     *  @return true, if lhs has less neighbors than rhs.
     */
    bool operator()(const MKVertex& lhs, const MKVertex& rhs) const
    {
        return lhs.neighbors < rhs.neighbors;
    }
};

/**
 *  Sorts the vertices of a subset of the vertices in a graph by their degree, descending.
 *
 *  @param G               Graph that contains the vertices and the vertex set.
 *  @param vertexSet       Vertex set to sort by degree.
 *  @param sortedVertexSet Array that will contain the vertices of the vertex set sorted by degree, descending.
 */
void VertexCoverLocalSearch::sortVerticesFromSetInGraphByDegree(vector<int> &vertexDegrees,
                                                                const Graph::vertex_set &vertexSet,
                                                                vector<Graph::vertex> &sortedVertexSet)
{
    /**
     *  We need to create a MKVertex object for each Graph::vertex to sort
     *  the vertices by degree. The vertices array will be the one that is
     *  filled with MKVertex elements and sorted afterwards.
     */
    vector<MKVertex> vertices = {};

    /**
     *  Create MKVertex objects for each vertex in the vertex set.
     */
    for (auto v : vertexSet)
    {
        MKVertex v_;
        v_.vertexID = v;
        v_.neighbors = vertexDegrees[v];
        vertices.push_back(v_);
    }

    /**
     *  Actually sort the vertices.
     */
    std::sort(vertices.begin(), vertices.end(), MKVertex());

    /**
     *  Write the sorted vertex IDs to the 'output' array.
     */
    for (int index = 0; index < vertices.size(); index++)
    {
        sortedVertexSet.push_back(vertices[index].vertexID);
    }
}

#pragma mark - Enumeration

void VertexCoverLocalSearch::addEdgesBetweenVerticesToGraph(vector<pair<Graph::vertex, Graph::vertex_set>> &vertices, MKGraph &G)
{
    for (auto p : vertices)
    {
        Graph::vertex v = p.first;
        Graph::vertex_set neighborsOfV = p.second;

        for (auto neighbor : neighborsOfV)
        {
            G.insert_edge(v, neighbor);
        }
    }
}

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
void VertexCoverLocalSearch::addBlackVertexToTreeGreedily(Graph::vertex &v,
                                                          const MKGraph &G,
                                                          const Graph::vertex_set &C,
                                                          Graph::vertex_set &T,
                                                          Graph::vertex_set &B,
                                                          Graph::vertex_set &TB,
                                                          Graph::vertex_set &P,
                                                          map<Graph::vertex, MKVertexDescendant> &descendants,
                                                          bool applyPruningRule2)
{
    /**
     *  We add v as a black vertex to T.
     */
    T.insert(v);

    /**
     *  v is a black vertex in the tree.
     */
    TB.insert(v);

    /**
     *  v must not be processed further.
     */
    B.insert(v);

    /**
     *  v is not a potential vertex anymore.
     */
    P.erase(v);

    /**
     *  Also all black neighbors of v must not be added to T ever, so we add them
     *  to B
     */
    Graph::vertex_set blackNeighborsOfV = G.neighbors(v) * C;

    /**
     *  Do not process them again.
     */
    B += blackNeighborsOfV;

    /**
     *  These black neighbors cannot be added to the tree anymore and therefore
     *  aren't potential vertices.
     */
    P -= blackNeighborsOfV;

    /**
     *  We only update the descendant structure if we're making
     *  use of it later.
     */
    if (applyPruningRule2)
    {
        /**
         *  We just added black vertex v to the set. All the
         *  ancestors of v now have one more black descendant.
         */
        MKVertexDescendant vDescendant;
        vDescendant.vertex = v;

        /**
         *  We're adding v greedily to T. We can add it as a
         *  child of any of its white neighbors that are in T.
         */
        Graph::vertex_set neighborsOfV = G.neighbors(v);
        Graph::vertex parentOfV = *neighborsOfV.begin();
        vDescendant.parent = parentOfV;
        descendants[v] = vDescendant;

        /**
         *  All ancestors of newly added vertex w have now one more
         *  black descendant.
         */
        descendants[parentOfV].addBlackDescendant(descendants);
    }
}

Graph::vertex_set VertexCoverLocalSearch::blackVerticesFromSet(const Graph::vertex_set &C,
                                                               const std::vector<bool> &vertexColors)
{
    Graph::vertex_set blackVertices = {};

    for (auto v : C)
    {
        if (vertexColors[v])
        {
            blackVertices.insert(v);
        }
    }

    return blackVertices;
}

Graph::vertex_set VertexCoverLocalSearch::whiteVerticesFromSet(const Graph::vertex_set &C,
                                                               const std::vector<bool> &vertexColors)
{
    Graph::vertex_set whiteVertices = {};

    for (auto v : C)
    {
        if (!vertexColors[v])
        {
            whiteVertices.insert(v);
        }
    }

    return whiteVertices;
}

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
bool VertexCoverLocalSearch::enumerateKExchangeSets(const MKGraph &G,
                                                    const Graph::vertex_set &C,
                                                    const std::vector<bool> &vertexColors,
                                                    int k,
                                                    Graph::vertex u,
                                                    Graph::vertex_set &T,   // Set of vertices in the tree
                                                    Graph::vertex_set &B,   // Set of black vertices that must not be processed further
                                                    Graph::vertex_set &W,   // Set of white vertices that must not be processed further
                                                    Graph::vertex_set &TB,  // Set of black vertices in the tree
                                                    Graph::vertex_set &TW,  // Set of white vertices in the tree
                                                    Graph::vertex_set &P,   // Set of black vertices that can still be added
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
                                                    long long &trianglesFoundWhileApplyingPruningRule3)
{
#ifdef DEBUG
    Debug("u = " << u << " ++++++++++++++++++++++++++++++++++++++++" << endl);
    Debug("[EnumParams]: #black = " << TB.size() << "\t#white = " << TW.size() << endl);

    Debug("[EnumParams]: T = ");
    MKGraph::printVertexSet(T);

    Debug("[EnumParams]: B = ");
    MKGraph::printVertexSet(B);

    Debug("[EnumParams]: W = ");
    MKGraph::printVertexSet(W);

    Debug("[EnumParams]: TB = ");
    MKGraph::printVertexSet(TB);

    Debug("[EnumParams]: P = ");
    MKGraph::printVertexSet(P);

    Debug("[EnumParams]: Descendants:" << endl);
    for (auto iterator = descendants.begin(); iterator != descendants.end(); iterator++)
    {
        Debug("[EnumParams]: \t" << iterator->first << ": #black: " << iterator->second.blackDescendants << ", #white: " << iterator->second.whiteDescendants << endl);
    }
#endif

    /**
     *  Whenever we enter this method, another node has been added to the search tree.
     */
    searchNodeCount += 1;

#pragma mark Cancelling Condition

    /**
     *  If these conditions are met, we found a valid exchange set.
     *
     *  Connected vertex set T is valid, if the number of black vertices
     *  is one more than the number of white vertices, while their sum
     *  is k.
     */
    if (T.size() == k &&
        TB.size() == (k / 2) + 1 &&
        TW.size() == (k / 2))
    {
        Debug("Found a valid exchange set!" << endl);
        return true;
    }
    /**
     *  We cannot find a valid vertex set T, if there are too many black or
     *  white vertices.
     */
    else if (TB.size() > (k / 2) + 1 ||
             TW.size() > (k / 2))
    {
        return false;
    }

    /**
     *  If the number of white vertices in T is k/2, we cannot add any more of them.
     *  Therefore, all black vertices that still have white neighbors that are not in T,
     *  can be added to the set B (not to be processed again), as adding them to T would
     *  require to add more white vertices to T as well. Additionally, these black
     *  vertices are removed from the set of potential black vertices P.
     */
    if (TW.size() == (k / 2))
    {
        Graph::vertex_set verticesThatCanNotBeAdded = {};

        /**
         *  Iterate the potential black vertices that could still be added.
         */
        for (auto v : P)
        {
            /**
             *  We can only check vertices that are in G. Here we're iterating
             *  the vertices in C. Some vertices might have been removed
             *  from the graph earlier!
             *  Additionally, we cannot add any black vertex that already has
             *  black neighbors in T.
             */
            if (G.find(v) != G.end())
            {
                /**
                 *  Here we check whether v has black neighbors in T, or white
                 *  neighbors that are not in T. In both cases v cannot be
                 *  added to T.
                 */
                Graph::vertex_set neighborsOfV = G.neighbors(v);

                Graph::vertex_set whiteNeighborsOfVThatAreNotInT = whiteVerticesFromSet(neighborsOfV, vertexColors);
                whiteNeighborsOfVThatAreNotInT -= TW;

                Graph::vertex_set blackNeighborsOfVThatAreInT = blackVerticesFromSet(neighborsOfV, vertexColors) * TB;

                if (whiteNeighborsOfVThatAreNotInT.size() > 0 ||
                    blackNeighborsOfVThatAreInT.size() > 0)
                {
                    verticesThatCanNotBeAdded.insert(v);
                }
            }
        }

        /**
         *  The vertices that cannot be added to T are now added to B (set of
         *  vertices that must not be processed again.) and removed from P (the
         *  set of potential black vertices.)
         */
        B += verticesThatCanNotBeAdded;
        P -= verticesThatCanNotBeAdded;
    }

    /**
     *  NOTE: The third cancelling condition (checking whether the black vertices
     *  in T are an independent set) is checked for implicitely, before
     *  entering another recursive step, as that will be faster.
     */

    /**
     *  First we get the neighbors of the current pivot vertex u.
     */
    Graph::vertex_set N_u = G.neighbors(u);

    /**
     *  Remove the vertices in B from the neighbors of pivot u, as they
     *  are not meant to be processed anymore.
     */
    N_u = N_u * P;

#pragma mark Pruning Rule 1

    /**
     *  We're no applying pruning rule 1, if enabled.
     */
    if (applyPruningRule1)
    {
        /**
         *  Iterate the black vertices that we're about to branch on and check
         *  whether some of the could be added greedily.
         */
        for (auto w : N_u)
        {
            /**
             *  We must not add more black vertices greedily, than are allowed.
             *  Therefore, we stop adding black vertices greedily if the
             *  threshhold is reached.
             */
            if (TB.size() == (k / 2) + 1)
            {
                break;
            }

            /**
             *  Initially, we assume to that the rule can be applied and check
             *  whether this still holds when checking for the conditions of
             *  the pruning rule.
             */
            bool rule1CanBeApplied = true;

            /**
             *  Now we're looking for a neighbor that does not met one of the
             *  criteria of the rule. If one is found, the rule cannot be applied.
             */
            Graph::vertex_set N_w = G.neighbors(w);

            for (auto x : N_w)
            {
                /**
                 *  If x is not white in T, one of the other two criteria
                 *  must be met.
                 */
                if (!includes_elm(TW, x))
                {
                    /**
                     *  Both of the remaining criteria require x to be in C
                     *  but not in T. Which means that x is a potential vertex.
                     */
                    if (!includes_elm(P, x))
                    {
                        rule1CanBeApplied = false;
                        break;
                    }
                    else // x is in C - T.
                    {
                        /**
                         *  Get all neighbors of neighbor x of w. (Except w of course)
                         */
                        Graph::vertex_set N_x = G.neighbors(x);
                        N_x.erase(w);

                        /**
                         *  If one of these neighbors is already a black vertex
                         *  in T, then x will never be removed.
                         *
                         *  We will instead check whether no neighbor is black
                         *  in T
                         */
                        if ((TB * N_x).size() == 0) // x has no black neighbor in T
                        {
                            /**
                             *  If the last criterium is not met either,
                             *  the pruning rule cannot be applied.
                             *
                             *  The last criterium is that all neighbors
                             *  of x are already in B (the set of vertices
                             *  that must not be processed further.) These
                             *  vertices will not be processed and they
                             *  didn't add x into T. Therefore x will never
                             *  be added to T and must be considered when
                             *  making the decision for w.
                             *
                             *  We will instead check if |N_x setminus B| > 0
                             *  which means that there are neighbors of x
                             *  that are not already in B.
                             */
                            if ((N_x - B).size() > 0)
                            {
                                /**
                                 *  At this point none of the requirements
                                 *  of the pruning rule where met, which
                                 *  means that the pruning rule cannot
                                 *  be applied.
                                 */
                                rule1CanBeApplied = false;
                                break;
                            }
                        }
                    }
                }
            }

            /**
             *  At this point the requirements for all neighbors of w
             *  have been checked. If the rule an be applied, we don't have
             *  decide whether to add w or not, but can be sure to add it.
             */

            if (rule1CanBeApplied)
            {
                Info("Pruning Rule 1 has been applied." << endl);

                addBlackVertexToTreeGreedily(w,
                                             G,
                                             C,
                                             T,
                                             B,
                                             TB,
                                             P,
                                             descendants,
                                             applyPruningRule2);

                /**
                 *  Pruning rule 1 has just been applied.
                 */
                pruningRule1NumberOfApplications += 1;
            }
        }

        /**
         *  We might just have added a few vertices greedily to T. Those vertices
         *  must not be processed again when branching.
         */
        N_u = N_u * P;
    }

#pragma mark Pruning Rule 3

    if (applyPruningRule3)
    {
        if (TW.size() == (k / 2))
        {
            /**
             *  We just reached the maximum number of white vertices in the exchange
             *  set. We check whether those have enough potential black neighbors to
             *  contribute to the set.
             *
             *  The number of black neighbors (that are not in the tree yet) of all
             *  white vertices in the set must be at least k/2+1. Else we cannot
             *  reach the size of the exchange set. More precisely: Since these
             *  black vertices might not be an independent set, we determine the
             *  size of a greedy matching of these black neighbors. If the number
             *  of isolated black vertices + the size of that greedy matching
             *  cannot increase the number of black vertices in T to k/2+1, we
             *  can stop diving deeper into the recursion here.
             */

            /**
             *  Now we get the union of all black neighbors of the white vertices in T.
             *  The white vertices are an independent set. Therefore, their neighborhood
             *  consists of black vertices by definition.
             *
             *  Additionally, we make sure to only search potential vertices
             *  that could still be added to T.
             */
            Graph::vertex_set blackNeighborsOfAllWhiteVertices = G.neighbors(TW) * P;

            /**
             *  The number of black vertices that might still be added to the T.
             */
            unsigned int numberOfPotentialBlackVertices = 0;

            /**
             *  We now create a subgraph that contains the black vertices only
             *  and remove all isolated vertices. After wards we will find the
             *  size of a greedy matching in the remaining subgraph.
             */
            MKGraph greedyMatchingSubgraph = G.subgraph(blackNeighborsOfAllWhiteVertices);

            /**
             *  First we remove all vertices with degree 0, as they cannot be
             *  part of a greedy matching. These vertices can be added to T
             *  anyways.
             */
            Graph::vertex_set isolatedVertices = {};

            for (MKGraph::const_iterator p = greedyMatchingSubgraph.begin();
                 p != greedyMatchingSubgraph.end();
                 p++)
            {
                const MKGraph::vertex &v = MKGraph::node(p);

                if (greedyMatchingSubgraph.neighbors(v).size() == 0)
                {
                    isolatedVertices.insert(v);
                }
            }

            /**
             *  We want to know how often we found isolated vertices when
             *  applying pruning rule 3.
             */
            if (!isolatedVertices.empty())
            {
                isolatedVerticesFoundWhileApplyingPruningRule3 += 1;
            }

            /**
             *  All the isolated black vertices can be added to the tree T.
             */
            numberOfPotentialBlackVertices += isolatedVertices.size();

            /**
             *  If the number of black potential vertices, that are isolated from
             *  each other, is >= k/2 + 1, we don't even need to find a greedy
             *  matching as there are already enough black vertices that might
             *  be added. If that is not the case, we calculate the size of a
             *  greedy matching.
             */
            if (TB.size() + numberOfPotentialBlackVertices < (k / 2) + 1)
            {
                /**
                 *  Remove the isolated black vertices from the subgraph.
                 */
                greedyMatchingSubgraph.remove_vertex_set(isolatedVertices);

                /**
                 *  Greedily find triangles and edges. For each triangle we have
                 *  to add at least two vertices. For each edge we have to add
                 *  at least one vertex.
                 */

                int trianglesFound = 0;
                int maximumBlackVerticesThatCanBeAddedToTheTree = greedyMatchingSubgraph.sizeOfGreedyTriangleMatching(trianglesFound);

                if (trianglesFound > 0)
                {
                    trianglesFoundWhileApplyingPruningRule3 += 1;
                }

                /**
                 *  Check whether the isolated vertices and the greedy matching
                 *  together are enough to expand the current exchange set to a
                 *  valid one.
                 */
                if (TB.size() + numberOfPotentialBlackVertices + maximumBlackVerticesThatCanBeAddedToTheTree < (k / 2) + 1)
                {
                    /**
                     *  We cannot expand the current exchange set to a valid one, since
                     *  there are not enough black vertices that we can add now.
                     *
                     *  Therefore, we stop the recursion here.
                     */
                    pruningRule3NumberOfApplications += 1;
                    return false;
                }
            }

            if (!isolatedVertices.empty())
            {
                /**
                 *  At this point we couldn't apply pruning rule 3. However, if
                 *  there were isolated black vertices here, we can add them to T
                 *  greedily without branching. We can do this, since black vertices
                 *  that have white neighbors that are not in T, have been eliminated
                 *  before.
                 */
                for (auto v : isolatedVertices)
                {
                    /**
                     *  We later have to find the parent of v in the tree. This is
                     *  done by intersecting the neighbors of v with the white
                     *  vertices in the set. We know that no more white vertices
                     *  are added and therefore can add v to any white vertex we
                     *  see fit.
                     */
                    Graph::vertex_set neighborsOfV = G.neighbors(v);

                    /**
                     *  v is only in the list of isolated vertices, if it has a
                     *  white neighbor in T, therefore the check whether the
                     *  neighbors of v are an empty set is unnecessary. Better
                     *  safe than sorry though.
                     */
                    if (!neighborsOfV.empty() &&
                        (neighborsOfV * TB).empty() &&
                        (whiteVerticesFromSet(neighborsOfV, vertexColors) - TW).empty() && // (((neighborsOfV - C) - T).empty()) &&
                        TB.size() < (k / 2) + 1) // We must not add more black vertices than necessary!
                    {
                        addBlackVertexToTreeGreedily(v,
                                                     G,
                                                     C,
                                                     T,
                                                     B,
                                                     TB,
                                                     P,
                                                     descendants,
                                                     applyPruningRule2);
                    }
                }
            }
        }

        /**
         *  When we're done here, we might have added some vertices greedily and
         *  and therefore need to check whether the conditions are met, yet.
         */

        N_u = N_u * P;
    }

    /**
     *  Here we're checking the cancelling conditions again. After applying the
     *  pruning rules, some black vertices might have been added greedily, which
     *  means that the exchange set might now be valid.
     *
     *  If these conditions are met, we found a valid exchange set.
     *
     *  Connected vertex set T is valid, if the number of black vertices
     *  is one more than the number of white vertices, while their sum
     *  is k.
     */
    if (T.size() == k &&
        TB.size() == (k / 2) + 1 &&
        TW.size() == (k / 2))
    {
        Debug("Found a valid exchange set!" << endl);
        return true;
    }
    /**
     *  We cannot find a valid vertex set T, if there are too many black or
     *  white vertices.
     */
    else if (TB.size() > (k / 2) + 1 ||
             TW.size() > (k / 2))
    {
        return false;
    }

#pragma mark Branching

    /**
     *  Iterate the neighbors w of u that still need to be processed and recursively
     *  try to find a valid exchange set by adding w to T.
     */
    for (auto w : N_u)
    {
        /**
         *  The neighbors of black vertex w which will be added to the set.
         */
        Graph::vertex_set N_w = G.neighbors(w);

        /**
         *  In the next recursive step we add w as black vertex to T.
         *  Since the black vertices in T have to be an independent set,
         *  we check whether a neighbor of w is a black vertex in T alreday.
         *
         *  N_w cap TB is the set containing neighbors of w that
         *  are also black vertices in T. If this is set is not empty, we
         *  cannot add w to T. If it is, we make a recursive call.
         */
        if ((N_w * TB).empty())
        {
            /**
             *  w is a neighbor of u, but since we're just processeing u, it
             *  must not be processed again here.
             */
            N_w.erase(u);

            /**
             *  The neighborhood of black vertex w can only contain white vertices.
             *  Therefore we subtract all black vertices (vertex cover C)
             *  and all vertices that should not be processed again (B).
             *  Additionally we don't want to add any vertices again,
             *  that are already in the tree. Note that there might be white
             *  vertices, that are still pivot elements to be processed and are
             *  therefore not in B or C, yet. That is why we make sure we don't
             *  process vertices again, that are already in the tree.
             */
            N_w = whiteVerticesFromSet(N_w, vertexColors); // N_w -= C;
            N_w -= TW;

            /**
             *  We just processed black vertex w and don't want it to be processed
             *  again in later recursive steps.
             */
            B.insert(w);

            /**
             *  w is no longer a potential vertex from now on.
             */
            P.erase(w);

            /**
             *  Here we construct the parameters for the next recursive call.
             *  Since most of the sets we're now editing are only valid inside
             *  the recursive call, we copy the current sets and edit the copies.
             *  These sets are marked es "*prime".
             */

            /**
             *  Black vertex w is now added to T
             */
            Graph::vertex_set Tprime = T;
            Tprime.insert(w);

            /**
             *  The set of black vertices in T now contains w
             */
            Graph::vertex_set TBprime = TB;
            TBprime.insert(w);

            /**
             *  The white neighbors of w are now added to T as well.
             */
            Tprime += N_w;

            /**
             *  N_w is the white neighborhood of w. They were just added to T
             *  and are therefore white vertices in T.
             */
            Graph::vertex_set TWprime = TW;
            TWprime += N_w;

            /**
             *  The descendants data structure that is used to keep track of
             *  how many black/white decendants a vertex has in T.
             */
            map<Graph::vertex, MKVertexDescendant> descendantsPrime = descendants;

            /**
             *  We only update the descendant structure if we're making use of
             *  it later, that is, if pruning rule 2 should be applied.
             */
            if (applyPruningRule2)
            {
                /**
                 *  We just added black vertex w to the set. All the
                 *  ancestors of w now have one more black descendant.
                 */
                MKVertexDescendant wDescendant;
                wDescendant.vertex = w;

                /**
                 *  w is a black child of white pivot u.
                 */
                wDescendant.parent = u;
                descendantsPrime[w] = wDescendant;

                /**
                 *  White pivot u and all its ancestor now have one more black
                 *  descendant.
                 */
                descendantsPrime[u].addBlackDescendant(descendantsPrime);
            }

            /**
             *  Add the white neighbors of w (N_w) to the top of the stack,
             *  since they'll be our next pivot elements.
             *
             *  Those white vertices have also been added to the set as children
             *  of w. Therefore, the descendants map has to be updated.
             */
            stack<Graph::vertex> pivotStackPrime = pivotStack;
            for (auto wprime : N_w)
            {
                pivotStackPrime.push(wprime);

                /**
                 *  We only update the descendant structure if we're making use of
                 *  it later.
                 */
                if (applyPruningRule2)
                {
                    MKVertexDescendant wprimeDescendant;
                    wprimeDescendant.vertex = wprime;

                    /**
                     *  w is the black parent of wprime.
                     */
                    wprimeDescendant.parent = w;

                    descendantsPrime[wprime] = wprimeDescendant;

                    /**
                     *  Parent w and its ancestors now have one more white descendant.
                     */
                    descendantsPrime[w].addWhiteDescendant(descendantsPrime);
                }
            }

            /**
             *  Also all black neighbors of w must not be added to T ever, so we add them
             *  to Bprime. We didn't add them to B, because when returning from
             *  the recursion, these neighbors are available again.
             */
            Graph::vertex_set Bprime = B;
            Graph::vertex_set blackNeighborsOfW = blackVerticesFromSet(G.neighbors(w), vertexColors); // G.neighbors(w) * C;
            Bprime += blackNeighborsOfW;

            /**
             *  These black neighbors are not in the set of potential vertices
             *  anymore, as they have a black neighbor in the tree.
             */
            Graph::vertex_set Pprime = P - blackNeighborsOfW;

            /**
             *  Now we add the black potential neighbors of the added white vertices.
             */
            Graph::vertex_set newPotentialBlackVertices = G.neighbors(N_w);

            /**
             *  Only black vertices are potential.
             */
            newPotentialBlackVertices = blackVerticesFromSet(newPotentialBlackVertices, vertexColors); //newPotentialBlackVertices * C;
            newPotentialBlackVertices -= Bprime;

            Pprime += newPotentialBlackVertices;

            /**
             *  Check if there are pivot vertices left. If there are,
             *  take the top most and dive deeper into the recursion.
             */
            if (!pivotStackPrime.empty())
            {
                /**
                 *  Now we're taking the first vertex from the stack of
                 *  pivot elements. If black vertex w had white neighbors
                 *  that weren't in the tree yet, they will be our next
                 *  pivot vertices. (Depth first approach).
                 *
                 *  If w didn't have any white neighbours to contribute to
                 *  out set T, u will continue to be our pivot vertex.
                 */
                Graph::vertex uprime = pivotStackPrime.top();

                /**
                 *  The recursive call with adapted parameters.
                 */
                if (enumerateKExchangeSets(G,
                                           C,
                                           vertexColors,
                                           k,
                                           uprime,
                                           Tprime,
                                           Bprime,
                                           W,
                                           TBprime,
                                           TWprime,
                                           Pprime,
                                           pivotStackPrime,
                                           searchNodeCount,
                                           descendantsPrime,
                                           applyPruningRule1,
                                           pruningRule1NumberOfApplications,
                                           applyPruningRule2,
                                           pruningRule2NumberOfApplications,
                                           applyPruningRule3,
                                           pruningRule3NumberOfApplications,
                                           isolatedVerticesFoundWhileApplyingPruningRule3,
                                           trianglesFoundWhileApplyingPruningRule3))
                {
                    T = Tprime;
                    return true;
                }
            }
        }
#ifdef DEBUG
        else
        {
            Debug(u << " has a black neighbor in T. Not going deeper here." << endl);
        }
#endif
    }

    /**
     *  Here we're done processing pivot u. Therefore, u should be at the
     *  top of the pivot stack, since we've been processing it until now.
     *  Checking for emptiness of the pivot stack should not be necessary.
     *  However, we're better safe than sorry.
     */
    if (!pivotStack.empty())
    {
        /**
         *  If the current pivot element u and the next pivot uprime have
         *  different parents, then we're done enumerating all subtrees of the
         *  black parent of first pivot u. If pruning rule 2 should be applied
         *  we now have to check whether that black parent has positive
         *  contribution. If not, we do not continue diving deeper into the
         *  recursion.
         */
        Graph::vertex parentU = descendants[u].parent;

        /**
         *  At this point we have completely processed u. It is therefore
         *  removed from the pivot stack.
         */
        pivotStack.pop();

        /**
         *  We've completely processed the last pivot element. Therefore,
         *  we continue with the next pivot element, if there is one.
         *
         *  Note: the following lines describe the construction of the
         *  branches of the search tree where we decided not to add any
         *  black neighbor of white pivot w to the set T. This means that
         *  the structure of T will not change in the next recursive call
         *  and therefore the cancelling conditions will not be met, if they
         *  haven't been met until now. Essentially, this means that we're
         *  not missing a valid solution by not diving deeper into the
         *  recursion if there is no new pivot element to choose.
         */
        if (!pivotStack.empty())
        {
            /**
             *  This is the case that represents not adding any neighbor of the
             *  current pivot.
             *  Get the next white pivot element and enumerate recursively.
             */
            Graph::vertex uprime = pivotStack.top();

#pragma mark Pruning Rule 2

            /**
             *  Check whether pruning rule 2 should be applied.
             */
            if (applyPruningRule2)
            {
                Graph::vertex parentUprime = descendants[uprime].parent;

                /**
                 *  If the two consecutive white pivot vertices (u and uprime)
                 *  have different parents, then we just finished enumerating the
                 *  subtree of the black parent of u. We now check whether
                 *  that black parent has positive contribution. If not, we
                 *  don't step deeper into the recursion.
                 */
                if (parentU != parentUprime)
                {
                    MKVertexDescendant descendantsOfBlackParentOfU = descendants[parentU];

                    /**
                     *  If the black parent of u has >= white descendants than
                     *  black ones, then it does not have positive contribution.
                     *
                     *  Note that a vertex is not a descendant of itself in our
                     *  MKVertexDescendant structure. Therefore we have to add
                     *  +1 to the number of black descendants to get the
                     *  contribution of the subtree rooted at u.
                     */
                    if (descendantsOfBlackParentOfU.whiteDescendants - (descendantsOfBlackParentOfU.blackDescendants + 1) >= 0)
                    {
                        /**
                         *  Pruning rule 2 has now been applied.
                         */
                        pruningRule2NumberOfApplications += 1;

#ifdef DEBUG
                        Debug("Applying pruning rule 2." << endl);
#endif

                        /**
                         *  Black parent of u does not have positive
                         *  contribution. Therefore, we don't step deeper into
                         *  the recursion.
                         */
                        return false;
                    }
                }
            }

            /**
             *  Make sure old pivot element u is not processed again.
             */
            Graph::vertex_set Wprime = W;
            W.insert(u);

            /**
             *  If we didn't return until now, we didn't find a valid set in a child
             *  node of the search tree. Therefore, we return whatever we'll find
             *  when continuing the recursion with the next pivot.
             */
            return enumerateKExchangeSets(G,
                                          C,
                                          vertexColors,
                                          k,
                                          uprime,
                                          T,
                                          B,
                                          Wprime,
                                          TB,
                                          TW,
                                          P,
                                          pivotStack,
                                          searchNodeCount,
                                          descendants,
                                          applyPruningRule1,
                                          pruningRule1NumberOfApplications,
                                          applyPruningRule2,
                                          pruningRule2NumberOfApplications,
                                          applyPruningRule3,
                                          pruningRule3NumberOfApplications,
                                          isolatedVerticesFoundWhileApplyingPruningRule3,
                                          trianglesFoundWhileApplyingPruningRule3);
        }
        else
        {
            /**
             *  Our stack of pivot elements is empty, which means there is
             *  nothing more to process, but we didn't meet the requirements
             *  yet. Return false.
             */
            return false;
        }
    }
    else
    {
        /**
         *  Our stack of pivot elements is empty, which means there is
         *  nothing more to process, but we didn't meet the requirements
         *  yet. Return false.
         */
        return false;
    }
}

/**
 *  Given a vertex Cover C, this method tries to find a vertex cover Cprime
 *  that has one vertex less than C by exchanging at most k vertices.
 *
 *  @param G                    Graph in which the vertex cover is to be improved.
 *  @param C                    Vertex cover to be improved.
 *  @param kMax                 Maximum number of exchanges to be done.
 *  @param Cprime               Improved solution, or an empty set if a better solution could not be found.
 *  @param applyPruningRule1    If set to true the first pruning rule will be applied. Default is false.
 *  @param applyPruningRule2    If set to true the second pruning rule (Positive contribution of subtree) will be applied. Default is false.
 *
 *  @return true, when a better vertex set Cprime was found. False, else.
 */
bool VertexCoverLocalSearch::improveVertexCoverLocallyByOne(MKGraph &G,
                                                            vector<int> &vertexDegrees,
                                                            Graph::vertex_set &C,
                                                            int kMax,
                                                            Graph::vertex_set &Cprime,
                                                            int elementToStartSearchAt,
                                                            int &elementToContinueSearchAt,
                                                            bool applyPruningRule1,
                                                            bool applyPruningRule2,
                                                            bool applyPruningRule3,
                                                            bool performSanityCheck)
{
    Info("Trying to improve vertex cover of size " << C.size() << " by exchanging <= " << kMax << " vertices." << endl);

    /**
     *  First we define the colors of the vertices.
     *
     *  White vertices are marked as false in the array of vertex colors.
     *  By default all vertices are white.
     */
    size_t numberOfVertices = G.num_vertices();
    std::vector<bool> vertexColors(numberOfVertices, false);

    /**
     *  All vertices that are in the cover are
     */
    for (auto v : C)
    {
        vertexColors[v] = true;
    }

    /**
     *  We will check for all odd values of k < kMax whether we can find a
     *  better solution by searching the k-exchange-neighborhood.
     */
    for (int k = 1; k <= kMax; k += 2)
    {
        /**
         *  To improve the performance of the algorithm we remove vertices from
         *  the graph, that we know will never be added during an enumeration.
         *  These will be vertices v that we already considered to be the first
         *  black vertex in T. If that case did not bring a valid exchange set,
         *  there won't be any exchange set that contains v.
         *
         *  We make a copy of G to not destroy the original graph.
         */
        //MKGraph Gprime(G);
        vector<pair<Graph::vertex, Graph::vertex_set>> removedVerticesAndNeighbors;

        Info("k = " << k << endl);

        /**
         *  We iterate the vertices in the given vertex set C.
         *  Each vertex v is chosen to be the first black vertex of T, that
         *  will be removed from C. Performance will be improved by starting
         *  with vertex with high degree to eliminate many edges at the
         *  beginning.
         */
        vector<Graph::vertex> coverSortedByDegree = {};
        VertexCoverLocalSearch::sortVerticesFromSetInGraphByDegree(vertexDegrees, C, coverSortedByDegree);

        /**
         *  We need an editable version of our cover to improve performance.
         */
        Graph::vertex_set vertexCover(coverSortedByDegree.begin(), coverSortedByDegree.end());

        /**
         *  Used to keep track of how many nodes in the search tree have been
         *  processed.
         */
        long long searchNodeCount = 0;

        /**
         *  Used to keep track of how often the second and third pruning rule has been
         *  applied.
         */
        long long pruningRule1NumberOfApplications = 0;
        long long pruningRule2NumberOfApplications = 0;
        long long pruningRule3NumberOfApplications = 0;

        /**
         *  Keep track of how often isolated vertices were found. (Not how many!)
         */
        long long isolatedVerticesFoundByApplyingPruningRule3 = 0;

        /**
         *  Keep track of how often triangles were found. (Not how many!)
         */
        long long trianglesFoundWhileApplyingPruningRule3 = 0;

        /**
         *  Here we start measuring how long we take to find an improvement or fail.
         */
        std::chrono::high_resolution_clock::time_point computationStartTime = std::chrono::high_resolution_clock::now();

        Debug("Iterating vertices in original vertex cover. v in C:" << endl);

        /**
         *  Specifying the start index allows us to perform cyclic local search.
         *  The search loop then works as follows:
         *  As long as possible (regarding the size), we extract the
         *  'elementToStartSearch'-th element from the array of cover vertices
         *  sorted by degree. We then try to find an improvement starting at that
         *  element. If none can be found, this vertex will never be part of a
         *  valid exchange set and is therefore removed from the graph and the cover.
         *  If the size of the cover is less than the index of that element,
         *  we start extracting vertices from the beginning of the array and
         *  check whether they are part of a valid set.
         *
         *  If elementToStartSearch = -1 this means we don't have an element to
         *  start our search at, which means we start at the beginning of the
         *  array, so the index is 0.
         */
        int index = 0;
        if (elementToStartSearchAt >= 0)
        {
            auto it = find(coverSortedByDegree.begin(), coverSortedByDegree.end(), elementToStartSearchAt);

            if (it != coverSortedByDegree.end())
            {
                index = (int)distance(coverSortedByDegree.begin(), it);
            }
        }

        /**
         *  Note that we remove those vertices v from the graph that were alreday
         *  processed as the initial black vertex in T. For those we checked
         *  every possible exchange set that wanted to remove v but none of them
         *  were valid exchange sets.
         *
         *  By removing them from the graph, we can save some computation time.
         *  This means that once there are less than (k/2 + 1) vertex cover nodes
         *  left in the graph, we cannot find a valid k-exchange set. Therefore,
         *  we can stop the iteration once we've removed |C| - (k/2) vertices.
         */
        while(coverSortedByDegree.size() >= max((k/2), 1))
        {
            /**
             *  v is the first black vertex in our exchange set.
             */
            Graph::vertex v = coverSortedByDegree[index];

            Debug("Processing start vertex v = " << v << endl);

            /**
             *  The following map will store how many black and white
             *  descendants a vertex will have in the exchange sets, which
             *  allows for pruning rules to be applied.
             */
            map<Graph::vertex, MKVertexDescendant> descendants;

            /**
             *  v is the first vertex, which does not have descendants yet.
             */
            MKVertexDescendant vDescendants;
            vDescendants.vertex = v;
            vDescendants.parent = v; // v is its own parent, which implies that this is the top of the tree.

            descendants[v] = vDescendants;

            /**
             *  We start with black vertex v and its white vertices N_w \ C
             */
            Graph::vertex_set T = {v};

            /**
             *  The set of black vertices in T now consists of {v}.
             */
            Graph::vertex_set TB = {v};

            /**
             *  Add the white neighbors of v to the tree.
             */
            Graph::vertex_set neighborsOfV = G.neighbors(v);//Gprime.neighbors(v);
            Graph::vertex_set whiteNeighborsOfV = whiteVerticesFromSet(neighborsOfV, vertexColors);

            T += whiteNeighborsOfV;

            /**
             *  The set of white vertices in T now consists of the white neighbors
             *  of v.
             */
            Graph::vertex_set TW = whiteNeighborsOfV;

            stack<Graph::vertex> pivotStack;

            /**
             *  Add all white neighbors of v to the stack as they will be
             *  pivot vertices. Additionally, we create MKVertexDescendant
             *  objects for each vertex to be able to track all the number of
             *  their descendants.
             */
            for (auto w : whiteNeighborsOfV)
            {
                pivotStack.push(w);

                MKVertexDescendant wDescandents;
                wDescandents.vertex = w;

                /**
                 *  White vertex w is added as a child of v.
                 */
                wDescandents.parent = v;
                descendants[w] = wDescandents;

                /**
                 *  v now has one more white child.
                 */
                descendants[v].addWhiteDescendant(descendants);
            }

            /**
             *  Check if there are pivot elements to process. If not,
             *  we found a black vertex v, whose neighbors are all black
             *  as well. This vertex didn't have to be in the cover anyway.
             *  Remove it from the cover and return the improved cover.
             */
            if (!pivotStack.empty())
            {
                /**
                 *  Choose the first pivot element u.
                 */
                Graph::vertex u = pivotStack.top();

                /**
                 *  B is the set of vertices that should not be processed later.
                 *  The black neighbors of our start vertex v must never be
                 *  added to the set, since no two black vertices in the set may
                 *  be adjacent.
                 */
                Graph::vertex_set B = blackVerticesFromSet(neighborsOfV, vertexColors); // N_v cap C

                /**
                 *  Vertex v must not be processed either.
                 */
                B.insert(v);

                /**
                 *  We only add potential vertices as we see them. At the moment
                 *  the only potential vertices are the black neighbors of all
                 *  white neighbors that are not in the tree yet. (v is the only
                 *  black vertex in the tree at the moment.)
                 */
                Graph::vertex_set P = blackVerticesFromSet(G.neighbors(whiteNeighborsOfV), vertexColors);
                P.erase(v);

                /**
                 *  Currently there are no white vertices that must not be processed further.
                 */
                Graph::vertex_set W = {};

                /**
                 *  Now we start to recusively enumerate the exchange sets.
                 */
                if (enumerateKExchangeSets(G, //Gprime
                                           vertexCover,
                                           vertexColors,
                                           k,
                                           u,
                                           T,
                                           B,
                                           W,
                                           TB,
                                           TW,
                                           P,
                                           pivotStack,
                                           searchNodeCount,
                                           descendants,
                                           applyPruningRule1,
                                           pruningRule1NumberOfApplications,
                                           applyPruningRule2,
                                           pruningRule2NumberOfApplications,
                                           applyPruningRule3,
                                           pruningRule3NumberOfApplications,
                                           isolatedVerticesFoundByApplyingPruningRule3,
                                           trianglesFoundWhileApplyingPruningRule3))
                {
                    /**
                     *  Measuring how long we took to find a valid exchange set.
                     */
                    std::chrono::high_resolution_clock::time_point computationEndTime = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(computationEndTime - computationStartTime).count();

                    /**
                     *  Printing some information about the process.
                     */
                    Info("Succeeded after " << duration / 1000.0 << " seconds, searching " << searchNodeCount << " nodes." << endl);
                    if (applyPruningRule1)
                    {
                        Info("(Pruning Rule 'Greedy Neighbors' was applied " << pruningRule1NumberOfApplications << " times.)" << endl);
                    }
                    if (applyPruningRule2)
                    {
                        Info("(Pruning Rule 'Positive Contribution' was applied " << pruningRule2NumberOfApplications << " times.)" << endl);
                    }
                    if (applyPruningRule3)
                    {
                        Info("(Pruning Rule 'Greedy Matching' was applied " << pruningRule3NumberOfApplications << " times.)" << endl);
                        Info("(Isolated vertices were found " << isolatedVerticesFoundByApplyingPruningRule3 << " times.)" << endl);
                        Info("(Triangles were found " << trianglesFoundWhileApplyingPruningRule3 << " times.)" << endl);
                    }
                    Reduced(k << "," << duration / 1000.0 << "," << searchNodeCount << "," << 1 << endl);

                    /**
                     *  T now contains the vertex set which defines which
                     *  vertices need to be exchanged.
                     */
#ifdef DEBUG
                    Debug("T = ");
                    MKGraph::printVertexSet(T);
#endif

                    /**
                     *  The black vertices are the ones to be removed from the
                     *  vertex cover. The white vertices are the ones to be added
                     *  to it.
                     */
                    Graph::vertex_set blackVertices = T * C;
                    Graph::vertex_set whiteVertices = T - blackVertices;

#ifdef DEBUG
                    Debug("Black vertices: ");
                    MKGraph::printVertexSet(blackVertices);
                    Debug("White vertices: ");
                    MKGraph::printVertexSet(whiteVertices);
#endif

                    /**
                     *  The new vertex cover is obtained by removing the black
                     *  vertices fromt the old cover and adding the white vertices.
                     */
                    Cprime = C - blackVertices + whiteVertices;

                    /**
                     *  Sanity check: test whether the newly created vertex
                     *  cover is indeed a vertex cover.
                     */
                    if (performSanityCheck)
                    {
                        if (G.isVertexCover(Cprime))
                        {
                            if (index + 1 < coverSortedByDegree.size())
                            {
                                elementToContinueSearchAt = coverSortedByDegree[index + 1];
                            }
                            else
                            {
                                elementToStartSearchAt = -1;
                            }
                            addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                            return true;
                        }
                        else
                        {
                            cout << "[Error]: Sanity check failed. The method found an alleged vertex cover that didn't cover all edges!" << endl;
                            addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                            return false;
                        }
                    }
                    else
                    {
                        addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                        return true;
                    }
                }
            }
            else
            {
                /**
                 *  The black vertex v that we wanted to remove from the
                 *  cover didn't have any neighbors that are not in the
                 *  cover already. Therefore v didn't have to be in the
                 *  cover in the first place. Remove v.
                 */
                Cprime = C;
                Cprime.erase(v);

                /**
                 *  Measuring how long we took to find a valid exchange set.
                 */
                std::chrono::high_resolution_clock::time_point computationEndTime = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(computationEndTime - computationStartTime).count();

                /**
                 *  Printing some information about the process.
                 */
                Info("Succeeded after " << duration / 1000.0 << " seconds, searching " << searchNodeCount << " nodes." << endl);
                if (applyPruningRule1)
                {
                    Info("(Pruning Rule 'Greedy Neighbors' was applied " << pruningRule1NumberOfApplications << " times.)" << endl);
                }
                if (applyPruningRule2)
                {
                    Info("(Pruning Rule 'Positive Contribution' was applied " << pruningRule2NumberOfApplications << " times.)" << endl);
                }
                if (applyPruningRule3)
                {
                    Info("(Pruning Rule 'Greedy Matching' was applied " << pruningRule3NumberOfApplications << " times.)" << endl);
                    Info("(Isolated vertices were found " << isolatedVerticesFoundByApplyingPruningRule3 << " times.)" << endl);
                    Info("(Triangles were found " << trianglesFoundWhileApplyingPruningRule3 << " times.)" << endl);
                }
                Reduced(k << "," << duration / 1000.0 << "," << searchNodeCount << "," << 1 << endl);

                /**
                 *  Sanity check: test whether the newly created vertex
                 *  cover is indeed a vertex cover.
                 */
                if (performSanityCheck)
                {
                    if (G.isVertexCover(Cprime))
                    {
                        if (index + 1 < coverSortedByDegree.size())
                        {
                            elementToContinueSearchAt = coverSortedByDegree[index + 1];
                        }
                        else
                        {
                            elementToStartSearchAt = -1;
                        }
                        addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                        return true;
                    }
                    else
                    {
                        cout << "[Error]: Sanity check failed. The method found an alleged vertex cover that didn't cover all edges!" << endl;
                        addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                        return false;
                    }
                }
                else
                {
                    addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
                    return true;
                }
            }

            /**
             *  Add this point we're done processing v as being the first
             *  black vertex in T. Therefore, future iterations don't have
             *  to consider v as a black vertex again.
             */

            removedVerticesAndNeighbors.push_back(make_pair(v, G.neighbors(v)));

            G.remove_vertex(v); //Gprime.remove_vertex(v);
            vertexCover.erase(v);

            /**
             *  Here we update the array, removing the vertex that we just
             *  searched. Afterwards the index is not changed, which indicates
             *  that we're continuing our search with the next vertex that took
             *  its place.
             */
            coverSortedByDegree.erase(coverSortedByDegree.begin() + index);

            if (coverSortedByDegree.size() <= index)
            {
                index = 0;
            }
        }

        /**
         *  Measuring the time it took we took to realize that we coudn't improve
         *  the cover here.
         */
        std::chrono::high_resolution_clock::time_point computationEndTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(computationEndTime - computationStartTime).count();
        Info("Failed after " << duration / 1000.0 << " seconds, searching " << searchNodeCount << " nodes." << endl);
        if (applyPruningRule1)
        {
            Info("(Pruning Rule 'Greedy Neighbors' was applied " << pruningRule1NumberOfApplications << " times.)" << endl);
        }
        if (applyPruningRule2)
        {
            Info("(Pruning Rule 'Positive Contribution' was applied " << pruningRule2NumberOfApplications << " times.)" << endl);
        }
        if (applyPruningRule3)
        {
            Info("(Pruning Rule 'Greedy Matching' was applied " << pruningRule3NumberOfApplications << " times.)" << endl);
            Info("(Isolated vertices were found " << isolatedVerticesFoundByApplyingPruningRule3 << " times.)" << endl);
            Info("(Triangles were found " << trianglesFoundWhileApplyingPruningRule3 << " times.)" << endl);
        }
        Reduced(k << "," << duration / 1000.0 << "," << searchNodeCount << "," << 0 << endl);

        addEdgesBetweenVerticesToGraph(removedVerticesAndNeighbors, G);
    }

    /**
     *  We haven't found a better vertex cover.
     */
    return false;
}

/**
 *  Calculates a vertex cover by iteratively removing the vertex with
 *  maximum degree from the graph.
 *
 *  @param cover A vertex cover of the graph.
 */
void VertexCoverLocalSearch::getGreedyVertexCoverFromGraph(const MKGraph &G, Graph::vertex_set &cover)
{
    MKGraph Gprime(G);

    while (Gprime.num_edges() > 0)
    {
        Graph::vertex maxDegV = Gprime.vertexWithMaxDegree();
        cover.insert(maxDegV);

        Gprime.remove_vertex(maxDegV);
    }
}

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
 *
 *  @return true, when a better vertex set Cprime was found. False, else.
 */
bool VertexCoverLocalSearch::improveVertexCoverLocally(const MKGraph &G,
                                                       Graph::vertex_set &C,
                                                       int kMax,
                                                       Graph::vertex_set &Cprime,
                                                       bool applyPruningRule1,
                                                       bool applyPruningRule2,
                                                       bool applyPruningRule3,
                                                       bool performSanityCheck)
{
    /**
     *  We try to find a better solution than C by repeatedly performing
     *  a local search on the currently best set..
     */
    Graph::vertex_set minimumSet = C;
    Graph::vertex_set improvedSet = {};

    /**
     *  Create an editable copy of the graph since "improveVertexCoverLocallyByOne"
     *  does not treat the graph as a constant.
     */
    MKGraph Gprime(G);

    /**
     *  Collect the degrees of the vertices in the graph. This information
     *  does not change but will be required very often so we don't want
     *  to generate it everytime.
     */
    vector<int> vertexDegrees = Gprime.vertexDegrees();

    /**
     *  We use this variable to describe at which element in the cover we want
     *  our next search to start. Initially the value is set to -1 which indicates
     *  that we start at index 0.
     *
     *  When an improvement is found the variable is overwritten with the element
     *  that would have been the next to search. This value is then passed as
     *  the next element to start the search at.
     */
    int elementToStartAt = -1;

    /**
     *  The reduced output will print the relevant information in a csv format.
     *  Here we're printing the headers.
     */
    Reduced("k,duration,searchNodes,improvementFound" << endl);

    while (VertexCoverLocalSearch::improveVertexCoverLocallyByOne(Gprime,
                                                                  vertexDegrees,
                                                                  minimumSet,
                                                                  kMax,
                                                                  improvedSet,
                                                                  elementToStartAt,
                                                                  elementToStartAt,
                                                                  applyPruningRule1,
                                                                  applyPruningRule2,
                                                                  applyPruningRule3,
                                                                  performSanityCheck))
    {
        minimumSet = improvedSet;
        improvedSet = {};

        Info("++++++++++++++++++++++++++++++++++++++++++++++++++" << endl);
        Info("Overall vertex cover improvement currently at: " << C.size() - minimumSet.size() << endl);
    }

    /**
     *  At this point we were not able to improve the last minimum set we
     *  had.
     */
    Cprime = minimumSet;

    /**
     *  If the vertex cover we found is smaller than the one we were given
     *  the algorithm was successful and we return true. Return false, if we
     *  didn't find a better solution.
     */
    if (Cprime.size() < C.size())
    {
        return true;
    }
    else
    {
        return false;
    }
}
