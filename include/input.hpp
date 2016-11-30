//
//  input.hpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 12.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#ifndef input_hpp
#define input_hpp

// C++ Includes
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include "ngraph.hpp"

using namespace NGraph;
using namespace std;

class Input
{
public:
    static int getGraphFromFile(string inputFileName, Graph &G);
    
    static int getVertexSetFromFile(string inputFileName, Graph::vertex_set &S);
};

#endif /* input_hpp */
