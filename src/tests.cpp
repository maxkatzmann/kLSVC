//
//  tests.cpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 15.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#include "tests.hpp"
#include "vertexcoverlocalsearch.hpp"

void Tests::getTestGraph(MKGraph &G)
{
    G = MKGraph();
    
    G.insert_edge(1, 2);
    G.insert_edge(2, 3);
    G.insert_edge(1, 3);
    G.insert_edge(1, 4);
}

bool Tests::testIsVertexCover()
{
    MKGraph G;
    getTestGraph(G);
    
    Graph::vertex one = 1;
    
    Graph::vertex_set nonCover = {one};
    
    if (G.isVertexCover(nonCover))
    {
        cout << "[Failure]: isVertexCover() confirmed a vertex cover that shouldn't be one." << endl;
        return false;
    }
    
    Graph::vertex three = 3;
    Graph::vertex_set cover = {one, three};
    if (!G.isVertexCover(cover))
    {
        cout << "[Failure]: isVertexCover() did not recognize a vertex cover that should be one." << endl;
        return false;
    }
    
    cout << "[Success]: Tested isVertexCover() successfully." << endl;
    return true;
}

bool Tests::testGetGreedyVertexCover()
{
    MKGraph G;
    getTestGraph(G);
    
    Graph::vertex_set cover = {};
    VertexCoverLocalSearch::getGreedyVertexCoverFromGraph(G, cover);
    
    if (!G.isVertexCover(cover))
    {
        cout << "[Failure]: getGreedyVertexCover(...) returned a vertex cover that was not really a vertex cover." << endl;
        return false;
    }
    
    cout << "[Success]: Tested getGreedyVertexCover(...) successfully." << endl;
    return true;
}

#pragma mark - Greedy Matching tests

bool Tests::testGreedyMatching()
{
    MKGraph G = MKGraph();
    
    G.insert_edge(1, 2);
    G.insert_edge(3, 4);
    G.insert_edge(5, 6);
    
    G.set_undirected();
    
    if (G.sizeOfGreedyMatching() != 3)
    {
        cout << "[Failure]: sizeOfGreedyMatching() didn't find a matching of the correct size." << endl;
        return false;
    }
    
    MKGraph G2 = MKGraph();
    
    /* A 4-clique has a matching of <= 2. But the any greedy matching should
     * have size 2.
     */
    G2.insert_edge(1, 2);
    G2.insert_edge(1, 3);
    G2.insert_edge(1, 4);
    G2.insert_edge(2, 3);
    G2.insert_edge(2, 4);
    G2.insert_edge(3, 4);
    
    G2.set_undirected();
    
    if (G2.sizeOfGreedyMatching() != 2)
    {
        cout << "[Failure]: sizeOfGreedyMatching() didn't find a matching of the correct size." << endl;
        return false;
    }
    
    cout << "[Success]: Tested sizeOfGreedyMatching() successfully." << endl;
    return true;
}

#pragma mark - Test Everything

bool Tests::testEverything()
{
    if (!testIsVertexCover())
    {
        return false;
    }
    
    if (!testGetGreedyVertexCover())
    {
        return false;
    }
    
    if (!testGreedyMatching())
    {
        return false;
    }
    
    cout << "[Success]: All tests succeeded." << endl;
    return true;
}