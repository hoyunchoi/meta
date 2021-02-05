#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <fstream>

#include "../library-Git/Networks.hpp"
#include "../library-Git/stringFormat.hpp"
#include "../library-Git/CSV.hpp"

struct Node_Meta : public WNode
{
    //* Member variables
    std::map<std::string, Size> m_numState;
    Size m_nodeSize;

    //* Generator
    Node_Meta(){}
    Node_Meta(const Size& t_index) : WNode(t_index){}
    Node_Meta(const Size& t_index, const std::vector<std::string>& t_states) : WNode(t_index){
        for (const auto& state : t_states){
            m_numState[state] = 0;
        }
    }

    void setSize(const Size& t_nodesize, const std::string& t_state = "S"){
        m_nodeSize = t_nodesize;
        m_numState[t_state] = t_nodesize;
    }
};

/*
    SIS Model simulation on Meta population network
    S + I -> I + I with rate SI_II
    I -> S with rate I_S
*/
//* Simulate SIS model on meata population network in deterministic way
namespace SIS_meta{
    //* Network parameters
    std::string networkType;
    Size networkSize;
    Size meanDegree;
    unsigned long long population;

    //* SIS rate parameter
    const std::vector<std::string> t_states = {"S", "I"};
    double SI_II;
    double I_S;
    double diffusion;

    //* Parameters used for RK4
    Size seedSize;
    Size maxIteration;
    double deltaT;
    unsigned long long numS, numI;
    std::vector<Node_Meta> nodes;
    std::set<Size> reactingIndex;

    //* Initialize
    void initialize(const WNetwork& t_wnetwork, const std::vector<double>& t_rates){
        //! Set Network
        networkType = t_wnetwork.m_type;
        networkSize = t_wnetwork.m_size;
        meanDegree = t_wnetwork.m_meanDegree;
        population = 0;

        //! Initialize model
        nodes.clear(); nodes.reserve(networkSize);
        reactingIndex.clear();
        for (Size index=0; index<networkSize; ++index){
            Node_Meta node(index, t_states);
            node.m_neighbors = t_wnetwork.m_wadjacency[index];      //* Weighted link and weight of node
            population += (Size)node.m_neighbors.at(index);         //* Add weight of node itself to population
            node.m_neighbors.erase(index);                          //* Delete weight of node itself in neighbor
            node.setSize(t_wnetwork.m_wadjacency[index].at(index), "S");
            nodes.emplace_back(node);
        }

        //! Set rates
        SI_II = t_rates[0];
        I_S = t_rates[1];
        diffusion = t_rates[2];

        //! Seed the model
        Size seedNum = 0;
        for (int index=0; index<networkSize; ++index){
            //* If chosen seed is isolated, choose again
            if (nodes[index].m_nodeSize != 0){
                nodes[index].m_numState["I"] += 1;
                nodes[index].m_numState["S"] -= 1;
                reactingIndex.emplace(index);
                ++seedNum;
            }
            if (seedNum >= seedSize){
                break;
            }
        }
        numI = seedNum;
        numS = population - seedNum;
    }//* End of function SIS_meta::initialize

    //* Update one step for every nodes
    void syncUpdate(){
        //* Intra-Node process
        for (const Size& index : reactingIndex){
            const double infecProb = 1.0 - std::exp(-1.0 * SI_II * deltaT * nodes[index].m_numState["I"] / nodes[index].m_nodeSize);
            const Size newInfected = std::floor(nodes[index].m_numState["S"] * infecProb + 0.5);

            const double recoverProb = I_S * deltaT;
            const Size newRecovered = std::floor(nodes[index].m_numState["I"] * recoverProb + 0.5);

            nodes[index].m_numState["S"] += newRecovered - newInfected;
            nodes[index].m_numState["I"] += newInfected - newRecovered;
            numS += newRecovered - newInfected;
        }
        numI = population - numS;

        //* Inter-Node process
        std::vector<int> deltaNumS(networkSize, 0);
        std::vector<int> deltaNumI(networkSize, 0);
        for (const Size& index : reactingIndex){
            for (auto it = nodes[index].m_neighbors.begin(); it != nodes[index].m_neighbors.end(); ++it){
                const double exchangeRatio = diffusion * it->second * deltaT;
                const int exchangeS = std::floor(exchangeRatio * nodes[index].m_numState["S"] + 0.5);
                const int exchangeI = std::floor(exchangeRatio * nodes[index].m_numState["I"] + 0.5);

                deltaNumS[index] -= exchangeS; deltaNumS[it->first] += exchangeS;
                deltaNumI[index] -= exchangeI; deltaNumI[it->first] += exchangeI;
            }
        }
        for (Size index=0; index<networkSize; ++index){
            nodes[index].m_numState["S"] += deltaNumS[index];
            nodes[index].m_numState["I"] += deltaNumI[index];
            nodes[index].m_nodeSize += deltaNumS[index] + deltaNumI[index];

            //* Update reacting Index
            if (nodes[index].m_numState["I"] > 0){
                reactingIndex.emplace(index);
            }
            else{
                reactingIndex.erase(index);
            }
        }
    }//* End of function SIS_meta::syncUpdate

    //* Run SIS_meta algorithm and save
    void syncRun(const WNetwork& t_wnetwork, const std::vector<double>& t_rates, const std::string& t_fileName = "", const std::string& t_totalFileName = "", const bool t_append = false){
        //* Open file and append
        std::ofstream file, totalFile;
        t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
        t_append ? totalFile.open(t_totalFileName, std::ios_base::app) : file.open(t_totalFileName);

        //* Initialize model
        initialize(t_wnetwork, t_rates);

        //* Print Initial values
        totalFile << numS << "," << numI << "\n";
        for(const Node_Meta& node : nodes){
            file << node.m_numState.at("S") << "," << node.m_numState.at("I") << ",";
        }
        file << "\n";

        //* Run sync update and print current values
        for (int iter=0; iter<maxIteration; ++iter){
            syncUpdate();
            totalFile << numS << "," << numI << "\n";
            for(const Node_Meta& node : nodes){
                file << node.m_numState.at("S") << "," << node.m_numState.at("I") << ",";
            }
            file << "\n";
        }

        //* new line for new ensemble
        // file << "\n";
        // totalFile << "\n";

        //* Close file
        file.close();
        totalFile.close();
    }//* End of function SIS_meta::syncRun
}//* End of namespace SIS_meta

/*
    SIR Model simulation on Meta population network
    S + I -> I + I with rate SI_II
    I -> R with rate I_R
*/
//* Simulate SIR model on meata population network in deterministic way
namespace SIR_meta{
    //* Network parameters
    std::string networkType;
    Size networkSize;
    Size meanDegree;
    unsigned long long population;

    //* SIR rate parameter
    const std::vector<std::string> t_states = {"S", "I", "R"};
    double SI_II;
    double I_R;
    double diffusion;

    //* Parameters used for RK4
    Size seedSize;
    Size maxIteration;
    double deltaT;
    unsigned long long numS, numI, numR;
    std::vector<Node_Meta> nodes;
    std::set<Size> reactingIndex;

    //* Initialize
    void initialize(const WNetwork& t_wnetwork, const std::vector<double>& t_rates){
        //! Set Network
        networkType = t_wnetwork.m_type;
        networkSize = t_wnetwork.m_size;
        meanDegree = t_wnetwork.m_meanDegree;
        population = 0;

        //! Initialize model
        nodes.clear(); nodes.reserve(networkSize);
        reactingIndex.clear();
        for (Size index=0; index<networkSize; ++index){
            Node_Meta node(index, t_states);
            node.m_neighbors = t_wnetwork.m_wadjacency[index];      //* Weighted link and weight of node
            population += (Size)node.m_neighbors.at(index);         //* Add weight of node itself to population
            node.m_neighbors.erase(index);                          //* Delete weight of node itself in neighbor
            node.setSize(t_wnetwork.m_wadjacency[index].at(index), "S");
            nodes.emplace_back(node);
        }

        //! Set rates
        SI_II = t_rates[0];
        I_R = t_rates[1];
        diffusion = t_rates[2];

        //* Seed the model
        Size seedNum = 0;
        for (int index=0; index<networkSize; ++index){
            //* If chosen seed is isolated, choose again
            if (nodes[index].m_nodeSize != 0){
                nodes[index].m_numState["I"] += 1;
                nodes[index].m_numState["S"] -= 1;
                reactingIndex.emplace(index);
                ++seedNum;
            }
            if (seedNum >= seedSize){
                break;
            }
        }
        numI = seedNum;
        numS = population - seedNum;
        numR = 0;
    }//* End of function SIR_meta::initialize

    //* Update one step for every nodes
    void syncUpdate(){
        //* Intra-Node process
        for (const Size& index : reactingIndex){
            const double infecProb = 1.0 - std::exp(-1.0 * SI_II * deltaT * nodes[index].m_numState["I"] / nodes[index].m_nodeSize);
            const Size newInfected = std::floor(nodes[index].m_numState["S"] * infecProb + 0.5);

            const double recoverProb = I_R * deltaT;
            const Size newRecovered = std::floor(nodes[index].m_numState["I"] * recoverProb + 0.5);

            nodes[index].m_numState["S"] -= newInfected;
            nodes[index].m_numState["I"] += newInfected - newRecovered;
            nodes[index].m_numState["R"] += newRecovered;
            numS -= newInfected;
            numR += newRecovered;
        }
        numI = population - numS - numR;

        //* Inter-Node process
        std::vector<int> deltaNumS(networkSize, 0);
        std::vector<int> deltaNumI(networkSize, 0);
        std::vector<int> deltaNumR(networkSize, 0);
        for (const Size& index : reactingIndex){
            for (auto it = nodes[index].m_neighbors.begin(); it != nodes[index].m_neighbors.end(); ++it){
                const double exchangeRatio = diffusion * it->second * deltaT;
                const int exchangeS = std::floor(exchangeRatio * nodes[index].m_numState["S"] + 0.5);
                const int exchangeI = std::floor(exchangeRatio * nodes[index].m_numState["I"] + 0.5);
                const int exchangeR = std::floor(exchangeRatio * nodes[index].m_numState["R"] + 0.5);

                deltaNumS[index] -= exchangeS; deltaNumS[it->first] += exchangeS;
                deltaNumI[index] -= exchangeI; deltaNumI[it->first] += exchangeI;
                deltaNumR[index] -= exchangeR; deltaNumR[it->first] += exchangeR;
            }
        }
        for (Size index=0; index<networkSize; ++index){
            nodes[index].m_numState["S"] += deltaNumS[index];
            nodes[index].m_numState["I"] += deltaNumI[index];
            nodes[index].m_numState["R"] += deltaNumR[index];
            nodes[index].m_nodeSize += deltaNumS[index] + deltaNumI[index] + deltaNumR[index];

            //* Update reacting Index
            if (nodes[index].m_numState["I"] > 0){
                reactingIndex.emplace(index);
            }
            else{
                reactingIndex.erase(index);
            }
        }
    }//* End of function SIR_meat::syncUpdate

    //* Run SIR_meta algorithm and save
    void syncRun(const WNetwork& t_wnetwork, const std::vector<double>& t_rates, const std::string& t_fileName = "", const std::string& t_totalFileName = "", const bool t_append = false){
        //* Open file and append
        std::ofstream file, totalFile;
        t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
        t_append ? totalFile.open(t_totalFileName, std::ios_base::app) : totalFile.open(t_totalFileName);

        //* Initialize model
        initialize(t_wnetwork, t_rates);

        //* Print Initial values
        totalFile << numS << "," << numI << "," << numR << "\n";
        for(const Node_Meta& node : nodes){
            file << node.m_numState.at("S") << "," << node.m_numState.at("I") << "," << node.m_numState.at("R") << ",";
        }
        file << "\n";

        //* Run sync update and print current values
        for (int iter=0; iter<maxIteration; ++iter){
            syncUpdate();
            totalFile << numS << "," << numI << "," << numR << "\n";
            for(const Node_Meta& node : nodes){
                file << node.m_numState.at("S") << "," << node.m_numState.at("I") << "," << node.m_numState.at("R") << ",";
            }
            file << "\n";
        }

        //* new line for new ensemble
        // file << "\n";
        // totalFile << "\n";

        //* Close file
        file.close();
        totalFile.close();
    }//* End of function SIR_meta::syncRun
}//* End of namespace SIR_meta

