//
//  main.cpp
//  kLSVC
//
//  Created by Maximilian Katzmann on 12.05.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#include <chrono>
#include <iostream>
#include <fstream>

#include "input.hpp"
#include "mkgraph.hpp"
#include "tests.hpp"
#include "vertexcoverlocalsearch.hpp"

using namespace NGraph;
using namespace std;

/**
 *  Prints information about the graph to standard output.
 *
 *  @param G The graph about which information should be printed.
 */
void printGraphInformation(MKGraph &G)
{
    Info("--------------------------------------------------" << endl);
    Info("Graph Information" << endl);
    Info("Number of nodes = " << G.num_vertices() << endl);
    Info("Number of edges = " << G.num_edges() << endl);
    Info("Max. Degree = " << G.maxDegree() << endl);
    Info("--------------------------------------------------" << endl);
//    Reduced("Nodes " << G.num_vertices() << endl);
//    Reduced("Edges " << G.num_edges() << endl);
//    Reduced("Max. Deg  " << G.maxDegree() << endl);
}

int main(int argc, const char * argv[])
{
    /**
     *  Options
     */
    string inputGraphFileName = "";
    string inputCoverFileName = "";
    string outputFileName = "";
    int kMax = 15;
    bool graphIsUndirected = true; // By default our graph is undirected.
    bool applyPruningRule1 = false;
    bool applyPruningRule2 = false;
    bool applyPruningRule3 = false;
    bool performSanityChecks = false;
    
    /**
     *  Parse command line arguments and update options
     */
    for (int index = 1; index < argc; index++)
    {
        if (string(argv[index]).compare("-graph") == 0)
        {
            if (index + 1 < argc)
            {
                inputGraphFileName = argv[index + 1];
            }
        }
        else if (string(argv[index]).compare("-directed") == 0)
        {
            graphIsUndirected = false;
        }
        else if (string(argv[index]).compare("-cover") == 0)
        {
            if (index + 1 < argc)
            {
                inputCoverFileName = argv[index + 1];
            }
        }
        else if (string(argv[index]).compare("-out") == 0)
        {
            if (index + 1 < argc)
            {
                outputFileName = argv[index + 1];
            }
        }
        else if (string(argv[index]).compare("-kMax") == 0)
        {
            if (index + 1 < argc)
            {
                kMax = stoi(argv[index + 1]);
            }
        }
        else if (string(argv[index]).compare("-pr1") == 0)
        {
            applyPruningRule1 = true;
        }
        else if (string(argv[index]).compare("-pr2") == 0)
        {
            applyPruningRule2 = true;
        }
        else if (string(argv[index]).compare("-pr3") == 0)
        {
            applyPruningRule3 = true;
        }
        else if (string(argv[index]).compare("-save") == 0)
        {
            performSanityChecks = true;
        }
        else if (string(argv[index]).compare("-help") == 0)
        {
            cout << "Usage: " << argv[0] << endl;
            cout << "\t-graph path/to/graphFile\t[Obligatory] Used to define the file containing the graph to be searched." << endl;
            cout << "\t-cover path/to/knownVertexCover\t[Optional] If set, the task is to find a better vertex cover. If not set, the task is to find some vertex cover." << endl;
            cout << "\t-out path/to/outputFile\t\t[Optional] If set, the result will be written to the specified file. If not set, the result will be printed to the console." << endl;
            cout << "\t-kMax int\t\t\t[Optional] An integer value > 0. If set, the exchange neighborhood of increasing size up until kMax is searched." << endl;
            cout << "\t-directed\t\t\t[Optional] Used to define whether the graph should be treated as directed or not. Default is undirected." << endl;
            cout << "\t-pr1\t\t\t\t[Optional] If set, the first pruning rule will be applied. False by default." << endl;
            cout << "\t-pr2\t\t\t\t[Optional] If set, the second pruning rule (Positive contribution of subtree) will be applied. False by default." << endl;
            cout << "\t-pr3\t\t\t\t[Optional] If set, the third pruning rule (White vertex max reached) will be applied. False by default." << endl;
            cout << "\t-help\t\t\t\tDisplays this help screen." << endl;
            cout << "\t-save\t\t\t\tEnables sanity checks for whether the found improved cover /really/ is a vertex cover. (decreases performance)" << endl;
            return 0;
        }
        else if (string(argv[index]).compare("-tests") == 0)
        {
            Info("--------------------------------------------------" << endl);
            Info("Performing tests:" << endl);
            Tests::testEverything();
            Info("--------------------------------------------------" << endl);
        }
    }
    
    if (inputGraphFileName == "")
    {
        cerr << "[Error]: Didn't provide input filename." << endl << "Usage: '" << argv[0] << " -graph /path/to/graphFile'" << endl << "Use: '" << argv[0] << " -help' for more information." << endl;
        return 1;
    }
    
    /**
     *  Creating the graph
     */
    MKGraph G;
    
    /**
     *  Setting whether the graph is undirected or not.
     */
    if (graphIsUndirected)
    {
        Info("Treating graph as undirected." << endl);
        G.set_undirected();
    }
    
    /**
     *  Reading the graph from the passed file.
     */
    Info("Reading graph from file: " << inputGraphFileName << endl);
    
    if (Input::getGraphFromFile(inputGraphFileName, G) != 0)
    {
        cout << "[Error]: The graph could not be read from the input file: " << inputGraphFileName << endl;
    }
    else
    {
        /**
         *  Printing information about the graph.
         */
        printGraphInformation(G);
        
        /**
         *  No input cover has been provided. Therefore we search one.
         */
        if (inputCoverFileName == "")
        {
            Graph::vertex_set greedyCover = {};
            
            VertexCoverLocalSearch::getGreedyVertexCoverFromGraph(G, greedyCover);
            
            if (G.isVertexCover(greedyCover))
            {
                Info("Found a greedy vertex cover with " << greedyCover.size() << " vertices!" << endl);
                
                /**
                 *  No output file has been provided. Therefore we print the result to the console.
                 */
                if (outputFileName == "")
                {
                    for (Graph::vertex_set::const_iterator v = greedyCover.begin(); v != greedyCover.end(); v++)
                    {
                        cout << *v << endl;
                    }
                }
                else
                {
                    // Write result to file.
                    fstream fileStream;
                    fileStream.open(outputFileName, fstream::out);
                    for (Graph::vertex_set::const_iterator v = greedyCover.begin(); v != greedyCover.end(); v++)
                    {
                        fileStream << *v << endl;
                    }
                    fileStream.close();
                    Info("Cover was written to output file: " << outputFileName << "." << endl);
                }
            }
            else
            {
                cout << "[Error]: Didn't find a greedy vertex cover!" << endl;
            }
        }
        else
        {
            /**
             *  The filename of a vertex cover was given. Read that cover into a
             *  vertex set.
             */
            Graph::vertex_set vertexCover = {};
            Input::getVertexSetFromFile(inputCoverFileName, vertexCover);
            /**
             *  Reading the graph from the passed file.
             */
            Info("Reading vertex cover from file: " << inputCoverFileName << endl);
            
            if (G.isVertexCover(vertexCover))
            {
                Info("Vertex cover from file is a valid cover." << endl);
                
                if (applyPruningRule1)
                {
                    Info("Pruning rule 1 will be applied whenever possible." << endl);
                }
                
                if (applyPruningRule2)
                {
                    Info("Pruning rule 2 (Positive Contribution of Subtree) will be applied whenever possible." << endl);
                }
                
                if (applyPruningRule3)
                {
                    Info("Pruning rule 3 (Greedy Matching) will be applied whenever possible." << endl);
                }
                
                if (performSanityChecks)
                {
                    Info("Sanity checks are enabled. Performance may be reduced." << endl);
                }
                
                /**
                 *  Now we try to find a better vertex cover by searching the neighborhood
                 *  of the given one.
                 */
                Graph::vertex_set betterVertexCover = {};
                std::chrono::high_resolution_clock::time_point computationStartTime = std::chrono::high_resolution_clock::now();
                if (VertexCoverLocalSearch::improveVertexCoverLocally(G,
                                                                      vertexCover,
                                                                      kMax,
                                                                      betterVertexCover,
                                                                      applyPruningRule1,
                                                                      applyPruningRule2,
                                                                      applyPruningRule3,
                                                                      performSanityChecks))
                {
                    std::chrono::high_resolution_clock::time_point computationEndTime = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(computationEndTime - computationStartTime).count();
                    Info("Found a better vertex cover in " << duration / 1000.0 << " seconds!" << endl);
                    Info("Size of the old cover: " << vertexCover.size() << ", size of the new cover: " << betterVertexCover.size() << ", improvement: " << vertexCover.size() - betterVertexCover.size() << endl);
                    
                    /**
                     *  No output file has been provided. Therefore we print the result to the console.
                     */
                    if (outputFileName == "")
                    {
                        for (Graph::vertex_set::const_iterator v = betterVertexCover.begin(); v != betterVertexCover.end(); v++)
                        {
                            cout << *v << endl;
                        }
                    }
                    else
                    {
                        // Write result to file.
                        fstream fileStream;
                        fileStream.open(outputFileName, fstream::out);
                        for (Graph::vertex_set::const_iterator v = betterVertexCover.begin(); v != betterVertexCover.end(); v++)
                        {
                            fileStream << *v << endl;
                        }
                        fileStream.close();
                        Info("Cover was written to output file: " << outputFileName << "." << endl);
                    }
                }
                else
                {
                    Info("Didn't find a better vertex cover." << endl);
                }
            }
            else
            {
                cout << "[Error]: Vertex cover from file is not a valid cover. An invalid cover cannot be improved!" << endl;
            }
        }
    }
    
    Info("Done." << endl);
    return 0;
}
