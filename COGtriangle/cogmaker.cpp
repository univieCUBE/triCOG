#include "cogmaker.h"
#include "iomanip"

#include <cassert>
#include <deque>
#include <queue>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>


// HitSet constructor
// Initialises HitSet class
HitSet::HitSet(	const double expect, 
				const double overlap,
				const std::map<int, std::vector<BlastHitData>>& blastHits,
				const std::map<int, std::string>& prot2org,
				const std::map<int, int>& prot2len,
				const std::map<int, Lse>& allLse,
				const std::map<int, int>&prot2lse,
				const std::string cName,
				const int startNum,
				const std::map<int, std::set<int>>& query2subject,
				std::shared_ptr<spdlog::logger> logger)
: organism("0")
, thresholds{overlap,expect}
, db{blastHits,prot2org,prot2len,allLse,prot2lse}
, filter{query2subject}
, logger(logger)
, cogName{cName}
, cogStartNum{startNum}
{   
    if (logger) {
        logger->info("HitSet initialized");
    } else {
		// Error if logger is null
        std::cerr << "Logger is zero" << std::endl;
    }
}

void HitSet::clear()
{
	organism = "0";
	lseHits.clear();
	conflict.clear();
}

void HitSet::insert(const std::set<int>& curLse) 
{
	// Check if set is empty
    if (curLse.empty()) {
        logger->error("Protein set is empty");
        throw std::invalid_argument("Cannot insert empty protein set");
    }
    // Log current LSE
    logger->debug("Processing LSE {}", *curLse.begin());
    // Warn if set is small
    if (curLse.size() < 3) {
        //logger->warn("insert(): only {} proteins (triangles need min 3)", curLse.size());
	}
		
	// Tests if first protein in curLse is in prot2org
	int id = *curLse.begin();
	auto p_s = db.prot2org.find(id);
	if (p_s == db.prot2org.end())
	{
		logger->error("Organism for protein {} not found", id);
    	throw std::out_of_range("Protein not in prot2org: " + std::to_string(id));
	}

	// Get the organism
	organism = p_s->second; // std::set organism for a given Lse

	// Find bidirectional best hits for each protein
	processLseProteins(curLse);

	// Decide whether to insert edge of LSE
	evaluateLseHits(curLse);

	// Go through conflict map and resolve by score
	resolveConflicts(curLse);
	
	// Log completion
    logger->debug("Finished processing {} proteins", curLse.size());
    logger->debug("Found {} best hits", lseHits.size());
	// Prepare for next LSE
	clear();
}

void HitSet::processLseProteins(const std::set<int>& curLse) {
	
	std::map<std::string, int> bet4q; // org, subject: BeTs for current query
	std::map<std::string, double> bestScore;

	for (const int& curQuery: curLse)
	{
		// Initialise iterator
		auto hitMapIter = db.blastHits.find(curQuery);
		// Check if query exists in blastHits
        if (hitMapIter == db.blastHits.end())
		{
			logger->warn("Protein {} has no BLAST hits", curQuery);
            continue; 
		} else 
		{
			bool use_filter = !filter.empty();
			auto filter_it = filter.end();
			// Apply filter logic
			if (use_filter) {
				filter_it = filter.find(curQuery);
				if (filter_it != filter.end()) {
					logger->debug("Query {} found in filter with {} allowed subjects", curQuery, filter_it->second.size());
				} else {
					logger->debug("Query {} NOT in filter - Skip this query", curQuery);
					continue;
				}
			}
			// Get vector of associated BLAST hits
			const std::vector<BlastHitData>& queryHits = hitMapIter->second;
			// Loop through each BLAST hit
			for (const BlastHitData& hit : queryHits)
			{
				int subject = hit.subject;

				// Validate subject exists
                auto prot2orgIt = db.prot2org.find(subject);
                auto prot2lseIt = db.prot2lse.find(subject);
                if (prot2orgIt == db.prot2org.end() ||
                    prot2lseIt == db.prot2lse.end())
                {
					 logger->warn("insert(): subject {} not found in database, skipping", subject);
					// Skip if subject not found
                    continue;  
                }
				// Apply custom filtering
				if (use_filter) { //&& filter_it != filter.end()
					const auto &allowed_subjects = filter_it->second;
					if (allowed_subjects.find(subject) == allowed_subjects.end()) {
						// Subject not in allowlist for this query
						continue;
					}
				}
				// Log exact hit
				logger->debug("Processing hit: query={}, subject={}, score={}, evalue={}, qStart={}, qEnd={}", 
				curQuery, hit.subject, hit.score, hit.evalue,hit.qStart, hit.qEnd);
				// Get organism, score and value
				std::string org = prot2orgIt->second; 
				double score = hit.score;  
				double evalue = hit.evalue; 

				// Skip self-hits and too high e-values
				if (org != organism && evalue <= thresholds.expect) // skip hits to self genome
				{
					// Parse alignment coordinates
					double qStart = hit.qStart; 
					double qEnd = hit.qEnd; 
					double sStart = hit.sStart; 
					double sEnd = hit.sEnd; 

					 // Validate coordinates are not negative
					if (qStart < 0 || qEnd < 0 || sStart < 0 || sEnd < 0) {
						logger->error("Negative coordinates for hit query={}, subject={}", curQuery, subject);
						logger->error("  qStart={}, qEnd={}, sStart={}, sEnd={}", qStart, qEnd, sStart, sEnd);
						throw std::domain_error("Coordinates cannot be negative");
					}
					
					// Validate start < end
					if (qStart > qEnd) {
						logger->error("Invalid query range qStart={} > qEnd={}", 
									qStart, qEnd);
						throw std::domain_error("Start coordinate > end coordinate");
					}
					if (sStart > sEnd) {
						logger->error("Invalid subject range sStart={} > sEnd={}", 
									sStart, sEnd);
						throw std::domain_error("Start coordinate > end coordinate");
					}

					// Calculate query and subject coverage and keep only if higher than overlap threshold
					if ((qEnd - qStart)/getLen(curQuery) >= thresholds.overlap &&
						(sEnd - sStart)/getLen(subject) >= thresholds.overlap)
					{
						// Check if we've seen hit to this genome before
						// If first hit: Store protein and score
						if (bet4q.find(org) == bet4q.end()) // first hit in this genome
						{
							bestScore.insert(std::pair<std::string, double>(org, score));
							bet4q.insert(std::pair<std::string, int>(org, subject));
						}
						// If not, compare score, keep better
						else if (score > (*bestScore.find(org)).second) // replace hit for better one
						{
							bestScore[org] = score;
							bet4q[org] = subject;
						}
					}						
				}	
			}

			// Iterate through best hits per genome
			for (const auto& entry : bet4q) {	
				int subject = entry.second;
				std::string org = entry.first;
				double score = (*bestScore.find(org)).second;
				// Find LSE for each protein
				int lse = (*db.prot2lse.find(subject)).second; // substitute subject protein with lse for this protein
				
				logger->debug("Added to lseHits: curQuery={}, subject={}, targetLse={}, org={}, score={}", 
							curQuery, subject, lse, org, score);

				// Insert Hsp object into multimap  keyed by LSE
				lseHits.insert(std::pair(lse, Hsp{lse,org,score}));
			}
			// Reset best-hit tracking for next protein in LSE
			bestScore.clear();
			bet4q.clear();
		}
	}
}

void HitSet::evaluateLseHits(const std::set<int>& curLse) {
	if (lseHits.empty()) return; // early return if no hits

	std::multimap<std::string, Hsp>::iterator p_last;
	int lseInserted = -1;
	double curLseSize = curLse.size();

	
	logger->debug("LSE {}: size = {} with {} hits", *curLse.begin(), curLse.size(), lseHits.size());
	// Iterate through all hits in current LSE
	for (const auto &hit : lseHits) {
		// Get LSE ID
		int lse = hit.first;

		// Print targets
		if (lse != lseInserted) {
        	logger->debug("   target LSE {} has count {}", lse, lseHits.count(lse));
		}
		
		// If target LSE was hit by more than half of the proteins
		if (lseHits.count(lse) > curLseSize/2 && lse != lseInserted)
		{
			// Store edge (query LSE ID, target LSE ID)
			sides.insert(std::pair<int, int>(*curLse.begin(), lse));
			// Track, to avoid duplication
			lseInserted = lse;
		}

		// Handle ambiguous cases (hits on exactly half)
		else if (lseHits.count(lse) == curLseSize/2)
		{
			std::string org = hit.second.subjectOrg;
			// If first into questionable
			if (lse != lseInserted)
			{
				// Store in conflict multimap
				p_last = conflict.insert(std::pair<std::string, Hsp>(org, hit.second));
				lseInserted = lse;
			}
			else // not first hit into questionable LSE
			{
				// Accumulate total score of 
				double score = hit.second.score;
				(*p_last).second.score += score; // compute total score for a questionable LSE
			}
		}
	} 
}

// Resolve conflict map
void HitSet::resolveConflicts(const std::set<int>& curLse) {
	// Early return
	if (conflict.empty()) return;

	// Go through conflict
	auto it = conflict.begin();
	while (it != conflict.end()) {
		const auto org = it->first;
		// Get all entries for this organism
		auto range = conflict.equal_range(org); // [first, second) for this org

		// If only one entry, add to sides
		// One hit with 50%, rest is spread out (no 50% -> not in this map)
		if (std::distance(range.first, range.second) == 1) {
			int lse = range.first->second.subjectLse;
			sides.insert({*curLse.begin(), lse});
		// If more than one, take the one with best score
		// This should be two (50 - 50)
		} else {
			auto first  = range.first;
			auto second = std::next(first);
			auto best   = (first->second.score < second->second.score) ? second : first;
			int lse = best->second.subjectLse;
			sides.insert({*curLse.begin(), lse});
		}
		// Jump to first entry of next organism
		it = range.second; 
	}
}

// Get organism in prot2org
std::string HitSet::getOrg(int prot)
{
	if (db.prot2org.find(prot) != db.prot2org.end())
	{
		return (*db.prot2org.find(prot)).second;
	}
	logger->warn("Protein {} not found in organism database", prot);
    return "-1";
}

// Find length in prot2len
int HitSet::getLen(int prot)
{
	if (db.prot2len.find(prot) != db.prot2len.end())
	{
		return (*db.prot2len.find(prot)).second;
	}
	logger->warn("Protein {} not found in length database", prot);
    return -1;
}


typedef std::set< std::pair<int, int> > edges_type;
typedef std::vector< std::vector<int> > vertex_map_type;
typedef std::vector<int>::const_iterator v_int_cit;

// Edge-based search: As each new edge is added, search for the new
// triangles that extend the COG that might result.

static inline void try_insert_edge(int v1, int v2, edges_type &cog_edges,
								   std::queue< std::pair<int, int> > &edge_queue)
{
	// Create pair of edges (forward, reverse)
	const std::pair<int, int> edge(v1, v2), edge_rev(v2, v1);

	// If edge not already in cog_edges
	if (cog_edges.find(edge) == cog_edges.end())
	{
		// Add edges to COG
		cog_edges.insert(edge);
		cog_edges.insert(edge_rev);
		// Add edge to BFS queue
		// Will later by processed by check_edge()
		edge_queue.push(edge);
	}
}

// Looks for all triangles containing one edge
static inline void check_edge(const std::pair<int, int> &edge,
							  const vertex_map_type &vertex_map,
							  const edges_type &symmetric_edges,
							  edges_type &cog_edges,
							  std::queue< std::pair<int, int> > &edge_queue)
{
	// Check: edge should already be in cog_edges
	assert(cog_edges.find(edge) != cog_edges.end());

	// Choose LSE with lower degree as seed
	// seed is edge endpoint having smaller degree (minor speedup)
	int seed = edge.first, nonseed = edge.second;
	if (vertex_map[edge.first].size() > vertex_map[edge.second].size())
		std::swap(seed, nonseed);

	// Iterate through every neighbour of seed
	for (v_int_cit p1=vertex_map[seed].begin();
		 p1 != vertex_map[seed].end(); p1++)
	{
		const int third = *p1;
		// If third vertex is first, continue
		if (third == nonseed)
			continue;
		// Check if triangle closes
		if (symmetric_edges.find(std::make_pair(third, nonseed))
			== symmetric_edges.end())
			continue;
		// Add new edges
		try_insert_edge(third, seed, cog_edges, edge_queue);
		try_insert_edge(third, nonseed, cog_edges, edge_queue);
	}
}

// An alternate (and perhaps better) strategy would be to interleave the
// second phase adds within the first phase, instead of explicitly switching
// between the two at coarse granularity, as is done here.

// Coordinates triangle finding
static inline void makeCOG(const int root_vertex, const int root_adj_vertex,
						   const vertex_map_type &vertex_map,
						   const edges_type &symmetric_edges,
						   edges_type &cog_edges, edges_type &processed)
{
	// Creates empty BFS queue
	std::queue< std::pair<int, int> > edge_queue;
	edge_queue.push(std::make_pair(root_vertex, root_adj_vertex));

	// In the first phase, try to form triangles from each edge in
	// edge_queue, adding new triangle edges to the std::queue.

	// Process edges until empty (FIFO order)
	while (not edge_queue.empty())
	  {
		// Remove oldest edge for processing
		const std::pair<int, int> edge = edge_queue.front();
		edge_queue.pop();

		// Mark edge as precessed
		processed.insert(edge);
		processed.insert(std::make_pair(edge.second, edge.first));

		// Check triangles from this edge
		check_edge(edge, vertex_map, symmetric_edges, cog_edges, edge_queue);
	  }
}

COGResults HitSet::makeCOGs()
{
    // Check if non-empty
    if (sides.empty()) {
        logger->warn("No edges found, no COGs to generate");
        COGResults results;
        return results;  // Return empty results
	}

	// Initialise results
	COGResults results;
	
	// Copies all edges from sides
	edges_type edges(sides.begin(), sides.end());

	// symmetric, non-self edges only
	edges_type symmetric_edges;

	// Get largest protein ID
	const int max_vertex = db.prot2org.rbegin()->first;

	// Checks for bi-directional edges (processes only one <)
	for (edges_type::iterator ep=edges.begin(); ep != edges.end(); ep++)
	{
		int from = ep->first, to = ep->second;
		if (edges.find(std::make_pair(to, from)) != edges.end()
			and from < to)
		{
			// Check if in prot2org and from different organisms
			std::map<int, std::string>::const_iterator from_p = db.prot2org.find(from);
			std::map<int, std::string>::const_iterator to_p = db.prot2org.find(to);
			if (from_p == db.prot2org.end()) {
                logger->error("Protein {} not found in prot2org", from);
                throw std::out_of_range("Edge protein not in database: " + std::to_string(from));
            }
            
            if (to_p == db.prot2org.end()) {
                logger->error("Protein {} not found in prot2org", to);
                throw std::out_of_range("Edge protein not in database: " + std::to_string(to));
            }
            
            if (from_p->second == to_p->second) {
                logger->error("Edge between proteins in same organism");
                logger->error("  from={} ({}), to={} ({})", from, from_p->second, to, to_p->second);
                throw std::logic_error("Should not have edges within same organism");
            }

			// Both directions in symmetric_edges
			symmetric_edges.insert(std::make_pair(from, to));
			symmetric_edges.insert(std::make_pair(to, from));
			
			// Save edges in results
			results.all_edges.push_back({from, to});
            results.all_edges.push_back({to, from});

			assert(0 <= from and from <= max_vertex);
			assert(0 <= to and to <= max_vertex);
		}
	}
	edges.clear();				// reclaim memory

	// vertex_map: std::map a vertex to a std::vector of adjacent vertices
	//             (same graph as symmetric_edges)
	vertex_map_type vertex_map(max_vertex+1);

	// Convert edge list into adjacency list
	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
		vertex_map[ep->first].push_back(ep->second);

	// True if an edge is already part of a COG or TWOG
	// (only contains ascending half-edge; v1 -> v2 where v1 < v2)
	edges_type processed;
	edges_type triprocessed;

	int cogN = cogStartNum;		// number assigned to cog in output

	// Iterate through all symmetric edges to find seed edges
	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
	{
		// Skip if reverse edge or edge already in processed
		const int root_vertex = ep->first, root_adj_vertex = ep->second;
		if (not (root_vertex < root_adj_vertex)
			or (processed.find(*ep) != processed.end()))
			continue;

		// Create empty COG and initialise seed edge in both directions
		edges_type cog_edges;
		// root_vertex--root_adj_vertex are endpoints of initial COG edge
		cog_edges.insert(std::make_pair(root_vertex, root_adj_vertex));
		cog_edges.insert(std::make_pair(root_adj_vertex, root_vertex));

		// Expand COG via triangle finding
		// now expand from that initial edge
		makeCOG(root_vertex, root_adj_vertex, vertex_map, symmetric_edges,
				cog_edges, processed);

		// Collect unique vertex IDs
		// if we found a COG, print it
		std::set<int> cog_vertices;	// vertices in cog_edges
		for (edges_type::const_iterator ep=cog_edges.begin();
			 ep != cog_edges.end(); ep++)
		  cog_vertices.insert(ep->first);

		// Filter for true COG
		// However, for true bidirectional graph 6 entries (two for each edge)
		if (cog_edges.size() >= 3)
		{
			// Write information to debug file
			for (edges_type::const_iterator ep=cog_edges.begin();
				 ep != cog_edges.end(); ep++)
			{
				triprocessed.insert(std::make_pair(ep->first, ep->second));
				triprocessed.insert(std::make_pair(ep->second, ep->first));
				int from = ep->first;
				int to = ep->second;
				// Create COG ID string
				std::ostringstream cogId;
				cogId << cogName << std::setw(5) << std::setfill('0') << cogN;
				
				// Store in results
				results.cog_edges.push_back({cogId.str(),from,to});
			}
			// Print final output, increment COG number
			returnCOG(cogN++, &cog_vertices, results);
		}
	}

#ifndef NO_TWOG_SINGLETON
	// Find bidirectional edges, that are not part of any triangle
	// print TWOGs (unfortunately a hack)
	for (edges_type::const_iterator ep=symmetric_edges.begin();
		 ep != symmetric_edges.end(); ep++)
	{
		const int root_vertex = ep->first, root_adj_vertex = ep->second;
		if (not (root_vertex < root_adj_vertex)
			or (triprocessed.find(*ep) != triprocessed.end()))
			continue;

		// Create minimal COG and return as COG
		edges_type cog_edges;
		// root_vertex--root_adj_vertex are endpoints of initial COG edge
		cog_edges.insert(std::make_pair(root_vertex, root_adj_vertex));
		cog_edges.insert(std::make_pair(root_adj_vertex, root_vertex));

		// If COG is found, add to results
		std::set<int> cog_vertices;	// vertices in cog_edges
		for (edges_type::const_iterator ep=cog_edges.begin();
			 ep != cog_edges.end(); ep++)
		  cog_vertices.insert(ep->first);
		
		returnCOG(cogN++, &cog_vertices, results);
	}

	// Collect all LSE IDs that were involved in any edge
	std::set<int> processed_vertices;
	for (edges_type::const_iterator p=processed.begin();
		 p != processed.end(); p++)
	{
		processed_vertices.insert(p->first);
		processed_vertices.insert(p->second);
	}

	// Find and output singletons
	for (std::map<int, Lse>::const_iterator p_lse = db.allLse.begin();
		 p_lse != db.allLse.end(); p_lse++)
	{
		const int vertex = p_lse->first;
		if (processed_vertices.find(vertex) == processed_vertices.end())
		{
			std::set<int> cog_vertices;
			cog_vertices.insert(vertex);
			
			returnCOG(0, &cog_vertices, results);
		}
	}
#endif
	logger->info("Completed COG generation: {} COG entries, {} edges in COGs", results.cog_entries.size(), results.cog_edges.size());
	return results;
}

// Save COG in results
void HitSet::returnCOG(int cogN, std::set<int>* cog, COGResults& results)
{
	std::string prot;
	std::set<int>::iterator p_lse;
	std::map<int, std::string>::iterator mp;

	for (std::set<int>::const_iterator p_cog = cog->begin();
		 p_cog != cog->end(); p_cog++)
	{
		const int lse = *p_cog;

		const std::set<int> lses = db.allLse.find(lse)->second.lse;
		for (std::set<int>::const_iterator p_lse = lses.begin();
			 p_lse != lses.end(); p_lse++)
		{
				// Create COG entry for this protein
				COGEntry entry;
				entry.protein = std::to_string(*p_lse);
				entry.organism = getOrg(*p_lse);
				entry.protein_name = "";
				entry.length = getLen(*p_lse);
				entry.start = 1;
				entry.end = getLen(*p_lse);
					if (cogN == 0)
				{
					entry.cog_id = std::nullopt;
				}
				else
				{
					std::ostringstream oss;
					oss << cogName << std::setw(5) << std::setfill('0') << cogN;
					entry.cog_id = oss.str();
				}
				results.cog_entries.push_back(entry);
		}
	}
}
