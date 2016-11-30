//
//  mkgraph.hpp
//  kLSVC
//
//  This class is used to add functionality to NGraph. This include method for
//  finding vertex covers of a graph, checking whether a given vertex set is a
//  vertex cover of a graph and the like.
//
//  Created by Maximilian Katzmann on 15.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#ifndef mkgraph_hpp
#define mkgraph_hpp

#ifdef DEBUG
#define Debug(x) cout << "[DEBUG]: " << x
#else
#define Debug(x)
#endif

#ifdef REDUCEDOUTPUT
#define Info(x)
#define Reduced(x) cout << x
#else
#define Reduced(x) 
#define Info(x) cout << "[Info]: " << x
#endif

#include "ngraph.hpp"
#include <stack>
#include <chrono>

using namespace NGraph;
using namespace std;

class MKGraph : public Graph
{
public:
    
    /**
     *  Prints a given vertex set to standard output.
     *
     *  @param vertexSet The vertex set to print.
     */
    static void printVertexSet(const Graph::vertex_set &vertexSet)
    {
        string setString = "{";
        for (Graph::vertex_set::const_iterator v = vertexSet.begin(); v != vertexSet.end(); v++)
        {
            setString.append(to_string(*v) + " ");
        }
        setString = setString.substr(0, setString.size() - 1);
        setString.append("}");
        
        cout << setString << endl;
    }
    
    /**
     *  Finds one vertex with the maximum degree. If there are multiple with the
     *  same maximum degree. One of the is returned.
     *
     *  @return A vertex with maximum degree.
     */
    vertex vertexWithMaxDegree() const
    {
        int maxDegree = 0;
        vertex vertexWithMaxDegree = 0;
        
        /**
         *  Iterate all vertices.
         */
        for (const_iterator p = begin(); p != end(); p++)
        {
            const vertex &v = node(p);
            
            /**
             *  Calculate the degree of the vertex.
             */
            unsigned int outDegree = out_degree(v);
            unsigned int inDegree = in_degree(v);
            
            unsigned int degree = outDegree + inDegree;
            
            /**
             *  Update maximum degree
             */
            if (degree > maxDegree)
            {
                maxDegree = degree;
                vertexWithMaxDegree = v;
            }
        }
        
        return vertexWithMaxDegree;
    }
    
    vector<int> vertexDegrees() const
    {
        vector<int> degrees((*this).num_nodes());
        
        /**
         *  Using the constructor to init the degrees vector with a given size
         *  raised an allocation error on linux systems. Therefore, we 'reserve'
         *  its size manually.
         */
        for (int index = 0; index < num_nodes(); index++)
        {
            degrees.push_back(0);
        }
        
        /**
         *  Iterate all vertices.
         */
        for (const_iterator p = begin(); p != end(); p++)
        {
            const vertex &v = node(p);
            
            /**
             *  Calculate the degree of the vertex.
             */
            unsigned int outDegree = out_degree(v);
            unsigned int inDegree = in_degree(v);
            
            unsigned int degree = outDegree + inDegree;
            
            degrees[v] = degree;
        }
        
        return degrees;
    }
    
    /**
     *  Determines the maximum degree of the graph.
     *
     *  @return The maximum degree of the graph.
     */
    unsigned int maxDegree() const
    {
        /**
         *  Get the vertex with maximum degree.
         */
        vertex maxDegV = vertexWithMaxDegree();
        return out_degree(maxDegV) + in_degree(maxDegV);
    }
    
    /**
     *  Finds all neighbors of a vertex in the graph.
     *
     *  @param v The vertex whose neighbors are to be obtained.
     *
     *  @return The neighbors of v.
     */
    const vertex_set neighbors(const vertex &v)const
    {
        if ((*this).is_directed())
        {
            /**
             *  If the graph is directed, only vertices that are adjaced to v
             *  with an outgoing edge from v are neighbors of v.
             */
            return (*this).out_neighbors(v);
        }
        else
        {
            /**
             *  If the graph is not directed all vertices that are adjacent with
             *  v are neighbors of v.
             */
            return (*this).out_neighbors(v) + (*this).in_neighbors(v);;
        }
    }
    
    /**
     *  Finds all neighbors of the vertices in the specified vertex set.
     *
     *  @param S Vertex set of the vertices whose neighbors are to be obtained.
     *
     *  @return A vertex set representing the union of the neighborhoods of all vertex in S.
     */
    const vertex_set neighbors(const vertex_set &S)const
    {
        vertex_set neighborOfS = {};
        for (auto v : S)
        {
            neighborOfS += (*this).neighbors(v);
        }
        
        return neighborOfS;
    }
    
    /**
     *  Checks whether the passed vertex set is a vertex cover.
     *
     *  @param S The supposed vertex cover.
     *
     *  @return true, if S is a vertex cover of the graph. false, if not.
     */
    bool isVertexCover(Graph::vertex_set &S) const
    {
        /**
         *  Copy the graph.
         */
        MKGraph G(*this);
        
        /**
         *  Remove the vertices of S and indident edges from G.
         */
        G.remove_vertex_set(S);
        
        /**
         *  Check if the graph without the vertex set does not have any edges.
         */
        return (G.num_edges() == 0);
    }
    
    /**
     *  Creates a graph containing only the vertices and edges between these
     *  vertices, that are in the passed vertex set.
     *
     *  @param A Set containing the vertices that should be in the subgraph.
     *
     *  @return A subgraph only conaining the vertices in the passed vertex set.
     */
    MKGraph subgraph(const vertex_set &A) const
    {
        MKGraph G;
        
        for (typename vertex_set::const_iterator p = A.begin(); p!=A.end(); p++)
        {
            const_iterator t = find(*p);
            if (t != end())
            {
                vertex_set new_in =  (A * in_neighbors(t));
                vertex_set new_out = (A * out_neighbors(t));
                
                G.insert_new_vertex_inout_list(*p, new_in, new_out);
            }
        }
        return G;
    }
    
    /**
     *  Calculates a greedy matching and returns its size.
     *
     *  @return The size of the greedy matching.
     */
    unsigned int sizeOfGreedyMatching()
    {
        /**
         *  Copy the graph, so we can remove vertices and edges without
         *  destroying the actual graph.
         */
        MKGraph G(*this);
        
        unsigned int size = 0;
        
        /**
         *  As long as there are edges in the graph...
         */
        while (G.num_edges() > 0)
        {
            /**
             *  ...we remove any edge and its endpoints, while increasing the
             *  size of the matching by one.
             *
             *  This is done by taking the first vertex with degree > 0 and
             *  removing it and one of its neighbors.
             */
            
            /**
             *  We remove all vertices with degree 0 as they will not be part of
             *  the matching anyway and we won't have to consider them twice.
             */
            Graph::vertex_set verticesWithDegreeZero = {};
            
            /**
             *  Iterate all vertices to find a vertex
             */
            for (const_iterator p = G.begin(); p != G.end(); p++)
            {
                const vertex &v = node(p);
                
                /**
                 *  Get the neighbors of the vertex.
                 */
                Graph::vertex_set neighbors = G.neighbors(v);
                
                if (!neighbors.empty())
                {
                    /**
                     *  We found a vertex that is endpoint of an edge. We will
                     *  consider this edge to be part of the matching and remove
                     *  both its endpoints from the graph.
                     *
                     *  Here we simply take the first neighbor we find.
                     */
                    vertex w = *neighbors.begin();
                    
                    /**
                     *  Remove both endporints.
                     */
                    G.remove_vertex(v);
                    G.remove_vertex(w);
                    
                    /**
                     *  We increased the matching by one.
                     */
                    size += 1;
                    
                    break;
                }
                else
                {
                    verticesWithDegreeZero.insert(v);
                }
            }
            
            /**
             *  Remove the vertices with zero degree so we don't have to process
             *  them again.
             */
            G.remove_vertex_set(verticesWithDegreeZero);
            
        }
        
        return size;
    }
    
    /**
     *  Tries to find triangles, greedily. If two vertices are not part of a 
     *  triangle they are removed and an edge is countet instead.
     *
     *  @return The size of the greedy matching.
     */
    
    /**
     *  Tries to find triangles, greedily. If two vertices are not part of a
     *  triangle they are removed and an edge is counted instead.
     *
     *  @param trianglesFound The number of triangles that were found.
     *
     *  @return Number of triangles and edges found.
     */
    unsigned int sizeOfGreedyTriangleMatching(int &trianglesFound)
    {
        /**
         *  Copy the graph, so we can remove vertices and edges without
         *  destroying the actual graph.
         */
        MKGraph G(*this);
        
        unsigned int size = 0;
        
        /**
         *  As long as there are edges in the graph...
         */
        while (G.num_edges() > 0)
        {
            /**
             *  ...we remove any triangle/edge and its endpoints, while increasing the
             *  size of the matching by one.
             *
             *  This is done by taking the first vertex with degree > 0 and
             *  removing it and one of its neighbors or two of its neighbors,
             *  if they are connected as well..
             */
            
            /**
             *  We remove all vertices with degree 0 as they will not be part of
             *  the matching anyway and we won't have to consider them twice.
             */
            Graph::vertex_set verticesWithDegreeZero = {};
            
            /**
             *  Iterate all vertices to find a vertex
             */
            for (const_iterator p = G.begin(); p != G.end(); p++)
            {
                const vertex &v = node(p);
                
                /**
                 *  Get the neighbors of the vertex.
                 */
                Graph::vertex_set neighbors = G.neighbors(v);
                
                if (!neighbors.empty())
                {
                    /**
                     *  We found a vertex that is endpoint of an edge. We will
                     *  consider this edge to be part of the matching and remove
                     *  both its endpoints from the graph.
                     *
                     *  Here we simply take the first neighbor we find.
                     */
                    vertex w = *neighbors.begin();
                    
                    Graph::vertex_set neighborsOfW = G.neighbors(w);
                    
                    Graph::vertex_set commonNeighbors = neighbors * neighborsOfW;
                    
                    /**
                     *  If v and w have a common neighbor, they form a triangle.
                     */
                    if (commonNeighbors.size() > 0)
                    {
                        /**
                         *  Remove one endpoint of the triangle. The other two
                         *  are removed later.
                         */
                        G.remove_vertex(*commonNeighbors.begin());
                        
                        trianglesFound += 1;
                    }
                    
                    /**
                     *  Remove both endpoints.
                     */
                    G.remove_vertex(v);
                    G.remove_vertex(w);
                    
                    /**
                     *  We increased the matching by one.
                     */
                    size += 1;
                    
                    break;
                }
                else
                {
                    verticesWithDegreeZero.insert(v);
                }
            }
            
            /**
             *  Remove the vertices with zero degree so we don't have to process
             *  them again.
             */
            G.remove_vertex_set(verticesWithDegreeZero);
            
        }
        
        return size;
    }
};

#endif /* mkgraph_hpp */
