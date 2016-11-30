//
//  input.cpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 12.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#include "input.hpp"
#include <stdexcept>

int Input::getGraphFromFile(string inputFileName, Graph &G)
{
    // Read content from file
    ifstream infile(inputFileName);

    // We keep track of the line number to be able to print where an error occured
    // when reading the file
    int lineNumber = 0;

    // Variable that will hold the content of each line of the input file.
    string line;

    // Read all lines in the file
    while (getline(infile, line))
    {
        lineNumber += 1;

        // Skip lines containging comments.
        if (line.size() < 1 || line.at(0) == '#' || line.at(0) == '%')
        {
            continue;
        }

        // Vector that will hold the labels of the nodes that are used in each line
        vector<string> nodeLabels;

        // Read the two labels from the line into the vector.
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(nodeLabels));

        // If a line contains more or less than 2 nodes, something is wrong. Print an error.
        if (nodeLabels.size() != 2)
        {
            ostringstream oss;
            oss << "[Error]: In line " << lineNumber << " the number of edges is != 2!";

            assert((oss.str(), nodeLabels.size() != 2));
        }
        else
        {
            // If no error occured, add the edge to the graph
            try {
                int node1 = stoi(nodeLabels[0]);
                int node2 = stoi(nodeLabels[1]);

                G.insert_edge(node1, node2);
            }
            catch (invalid_argument&)
            {
                cout << "[Error]: At least one of the nodes in line " << lineNumber << " was not an integer." << endl;
                return -1;
            }
        }
    }

    if (lineNumber == 0)
    {
        cout << "[Error]: The input file seemed to be empty!" << endl;
        return -1;
    }

    return 0;
}

int Input::getVertexSetFromFile(string inputFileName, Graph::vertex_set &S)
{
    // Read content from file
    ifstream infile(inputFileName);

    // We keep track of the line number to be able to print where an error occured
    // when reading the file
    int lineNumber = 0;

    // Variable that will hold the content of each line of the input file.
    string line;

    // Read all lines in the file
    while (getline(infile, line))
    {
        lineNumber += 1;

        // Skip lines containging comments.
        if (line.size() < 1 || line.at(0) == '#')
        {
            continue;
        }

        // If no error occured, add the edge to the graph
        try
        {
            int node = stoi(line);
            S.insert(node);
        }
        catch (invalid_argument&)
        {
            cout << "[Error]: The element in line " << lineNumber << " was not an integer." << endl;
            return -1;
        }
    }

    if (lineNumber == 0)
    {
        cout << "[Error]: The input file seemed to be empty!" << endl;
        return -1;
    }

    return 0;
}
