#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <memory>
#include <cstdlib>
#include <optional>
#include <spdlog/spdlog.h>




//#define NO_TWOG_SINGLETON

// Stores thresholds
struct Thresholds
{
	double	overlap; // minimum overlap between 2 sectors of the same COG to be joined
	double	expect; // e-value cut-off
};

// Stores LSE
struct Lse
{
	std::set<int> lse;
};

// Holds BLAST output
struct BlastHitData{
	int query;
	int subject;
	double qStart;
	double qEnd;
	double sStart;
	double sEnd;
	double evalue;
	double score;
};

// Holds input data
struct DB
{
	std::map<int, std::vector<BlastHitData>> blastHits;
	std::map<int, std::string>	prot2org;
	std::map<int, int>	prot2len;
	std::map<int, Lse> 	allLse; // every protein has lse, at least for itself (lse.size() = 1)
	std::map<int, int>	prot2lse;
};


struct Hsp
{
	int		subjectLse;
	std::string	subjectOrg;
	double	score;
};

// Stores one output protein
struct COGEntry
{
	std::string protein;
    std::string organism;
    std::string protein_name;
    int length;
    int start;
    int end;
    std::optional<std::string> cog_id;
};

// Holds all valid edges
struct AllEdges
{
	int from_id;
	int to_id;
};

// Holds all edges in cluster
struct COGEdges
{
	std::string cog_id;
	int from_id;
	int to_id;
};

// Holds result of programme
struct COGResults
{
	std::vector<COGEntry> cog_entries;
    std::vector<AllEdges> all_edges;      
    std::vector<COGEdges> cog_edges;      
};


class HitSet
{
friend class HitSetInspector; // for testing private members
public:
			HitSet	(const double expect,
					const double overlap,
					const std::map<int, std::vector<BlastHitData>>& blastHits,
					const std::map<int, std::string>&	prot2org,
					const std::map<int, int>& prot2len,
					const std::map<int, Lse>& allLse,
					const std::map<int, int>& prot2lse,
					const std::string cName,
					const int startNum,
					const std::map<int, std::set<int>>& query2subject,
					std::shared_ptr<spdlog::logger> logger = nullptr);

			void insert	(const std::set<int>& curLse);
			COGResults makeCOGs();
			void clear	();	
	
private:
	std::string	organism;
	Thresholds	thresholds;
	DB db;
	std::map<int, std::set<int>> filter;
	std::string	cogName;
	int	cogStartNum;
	std::shared_ptr<spdlog::logger> logger;
	std::multimap<int, Hsp>	lseHits; // for current LSE only: subjectLse, Hsp
	std::multimap<std::string, Hsp>	conflict; // subjectOrg, Hsp
	std::multimap<int, int>	sides; // apexQ, apexT
	std::set<int> processedLse;
	void processLseProteins(const std::set<int>& curLse);
	void evaluateLseHits(const std::set<int>& curLse);
	void resolveConflicts(const std::set<int>& curLse);
	std::string getOrg(int prot);
	int getLen(int prot);
	void returnCOG(int cogN, std::set<int>* cog, COGResults& results);
};
