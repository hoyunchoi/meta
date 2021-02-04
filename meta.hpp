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
    //* pre-defined parameters
    const std::string rootDirectory = "../data/epidemics/SIS_meta/";
    const std::vector<std::string> t_states = {"S", "I"};

    //* Network parameters
    std::string networkType;
    Size networkSize;
    Size meanDegree;

    //* SIS rate parameter
    double SI_II;
    double I_S;
    double diffusion;

    //* Parameters used for RK4
    Size seedSize;
    Size maxIteration;
    double deltaT;
    std::vector<Node_Meta> nodes;
    std::set<Size> reactingIndex;

    //* Initialize
    void initialize(const WNetwork& t_wnetwork, const std::vector<double>& t_rates){
        //! Set Network
        networkType = t_wnetwork.m_type;
        networkSize = t_wnetwork.m_size;
        meanDegree = t_wnetwork.m_meanDegree;

        //! Initialize model
        nodes.clear(); nodes.reserve(networkSize);
        reactingIndex.clear();
        for (Size index=0; index<networkSize; ++index){
            Node_Meta node(index, t_states);
            node.m_neighbors = t_wnetwork.m_wadjacency[index];      //* Weighted link and weight of node
            node.m_neighbors.erase(index);                          //* Delete weight of node itself in neighbor
            node.setSize(t_wnetwork.m_wadjacency[index].at(index), "S");
            nodes.emplace_back(node);
        }

        //! Set rates
        SI_II = t_rates[0];
        I_S = t_rates[1];
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
    }

    //* update one step for every nodes
    void syncUpdate(){
        //* Intra-Node process
        for (const Size& index : reactingIndex){
            const double infecProb = 1.0 - std::exp(-1.0 * SI_II * deltaT * nodes[index].m_numState["I"] / nodes[index].m_nodeSize);
            const Size newInfected = std::floor(nodes[index].m_numState["S"] * infecProb + 0.5);

            const double recoverProb = I_S * deltaT;
            const Size newRecovered = std::floor(nodes[index].m_numState["I"] * recoverProb + 0.5);

            nodes[index].m_numState["S"] += newRecovered - newInfected;
            nodes[index].m_numState["I"] += newInfected - newRecovered;
        }

        //* Inter-Node process
        std::vector<int> deltaNumI(networkSize, 0);
        for (const Size& index : reactingIndex){
            for (auto it = nodes[index].m_neighbors.begin(); it != nodes[index].m_neighbors.end(); ++it){
                const double exchangeRatio = diffusion * it->second * deltaT;
                const int exchangeI = std::floor(exchangeRatio * nodes[index].m_numState["I"] + 0.5);
                deltaNumI[index] -= exchangeI;
                deltaNumI[it->first] += exchangeI;
            }
        }
        reactingIndex.clear();
        for (Size index=0; index<networkSize; ++index){
            if (deltaNumI[index]){
                nodes[index].m_numState["I"] = std::min(nodes[index].m_nodeSize, std::max((Size)0, nodes[index].m_numState["I"] + deltaNumI[index]));
                nodes[index].m_numState["S"] = nodes[index].m_nodeSize - nodes[index].m_numState["I"];
            }
            //* Update reacting Index
            if (nodes[index].m_numState["I"] > 0){
                reactingIndex.emplace(index);
            }
        }
    }//* End of function SIS_meta::syncUpdate

    //* Run SIS_meta algorithm for ensemble size times
    void syncRun(const WNetwork& t_wnetwork, const std::vector<double>& t_rates, const std::string& t_fileName = ""){
        //* Open file and append
        std::ofstream file(t_fileName, std::ios_base::app);

        //* Initialize model
        initialize(t_wnetwork, t_rates);

        //* Print Initial values
        for(const Node_Meta& node : nodes){
            file << node.m_numState.at("S") << "\t" << node.m_numState.at("I") << "\t";
        }

        //* Run sync update and print current values
        for (int iter=0; iter<maxIteration; ++iter){
            syncUpdate();
            for(const Node_Meta& node : nodes){
                file << node.m_numState.at("S") << "\t" << node.m_numState.at("I") << "\t";
            }
        }

        //* new line for new ensemble
        file << "\n";
        file.close();
    }//* End of function SIS_meta::syncRun
}//* End of namespace SIS_meta