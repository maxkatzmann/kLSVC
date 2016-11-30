//
//  tests.hpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 15.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#ifndef tests_hpp
#define tests_hpp

#include <stdio.h>
#include "mkgraph.hpp"

using namespace std;

class Tests
{
private:
    static void getTestGraph(MKGraph &G);
    
public:
    
#pragma mark - VertexCover tests

    static bool testIsVertexCover();

    static bool testGetGreedyVertexCover();
    
#pragma mark - Greedy Matching tests
    
    static bool testGreedyMatching();
    
#pragma mark - Test Everything
    
    static bool testEverything();
};

#endif /* tests_hpp */
