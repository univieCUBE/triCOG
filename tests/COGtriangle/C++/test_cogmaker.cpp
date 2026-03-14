#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <spdlog/sinks/null_sink.h>
#include "cogmaker.h"

// Access to private members of HitSet for testing purposes
class HitSetInspector : public HitSet {
public:
    using HitSet::HitSet;
    
    const Thresholds& getThresholds() const {return thresholds; }
    std::string getCogName() const { return cogName; }
    int getCogStartNum() const { return cogStartNum; }
    const std::multimap<int, int>& getSides() const { return sides; }
};

TEST_CASE("HitSet Constructor: Validates initialization of all members", "[HitSet][init]") {
    
    // Get Mock data
    std::map<int, std::vector<BlastHitData>> mockBlastHits;
    std::map<int, std::string> mockProt2Org = {{1, "OrgA"}, {10, "OrgB"}};
    std::map<int, int> mockProt2Len = {{1, 150}, {10, 300}};
    std::map<int, Lse> mockAllLse;
    std::map<int, int> mockProt2Lse = {{1, 101}, {10, 202}};
    std::map<int, std::set<int>> mockQuery2Subject = {{1, {10}}};

    double expectVal = 0.001;
    double overlapVal = 0.6;
    std::string testName = "TEST_COG";
    int startNum = 500;

    // Initialise Hitset
    HitSetInspector inspector(
        expectVal, overlapVal, mockBlastHits, mockProt2Org, 
        mockProt2Len, mockAllLse, mockProt2Lse, testName, 
        startNum, mockQuery2Subject
    );

    // Check if members are initialized correctly
    SECTION("Thresholds are correctly mapped") {
        CHECK(inspector.getThresholds().expect == expectVal);
        CHECK(inspector.getThresholds().overlap == overlapVal);
    }

    SECTION("Cluster members are correctly assigned") {
        CHECK(inspector.getCogName() == "TEST_COG");
        CHECK(inspector.getCogStartNum() == 500);
    }
}

TEST_CASE("HitSet Constructor: Handles empty input gracefully", "[HitSet][init][edge_case]") {
    // Setup mock data
    std::map<int, std::vector<BlastHitData>> emptyHits;
    std::map<int, std::string> emptyP2O;
    std::map<int, int> emptyP2L;
    std::map<int, Lse> emptyLse;
    std::map<int, int> emptyP2Lse;
    std::map<int, std::set<int>> emptyQ2S;

    HitSetInspector inspector(
        1.0, 0.0, emptyHits, emptyP2O, emptyP2L, 
        emptyLse, emptyP2Lse, "", 0, emptyQ2S
    );

    // Check empty initialisation
    CHECK(inspector.getCogName().empty());
    CHECK(inspector.getCogStartNum() == 0);
}

TEST_CASE("HitSet::insert Basic insert check and filter logic", "[HitSet]") {
    // Setup mock data
    std::map<int, std::vector<BlastHitData>> blastHits;
    std::map<int, std::string> prot2org = {{100, "OrgA"}, {200, "OrgB"}};
    std::map<int, int> prot2len = {{100, 1000}, {200, 1000}};
    std::map<int, Lse> allLse = {{100, {{100}}}, {200, {{200}}}};
    std::map<int, int> prot2lse = {{100, 100}, {200, 200}};

    // One hit from 100 to 200
    blastHits[100] = {{100, 200, 0.0, 900.0, 0.0, 900.0, 1e-50, 500.0}};
    
    std::set<int> curLse = {100};
    double expect = 1e-10;
    double overlap = 0.5;
    auto test_logger = std::make_shared<spdlog::logger>("test", std::make_shared<spdlog::sinks::null_sink_st>());


    SECTION("Empty filter allows all hits") {
        std::map<int, std::set<int>> emptyFilter;
        HitSetInspector inspector(expect, overlap, blastHits, prot2org, prot2len, allLse, prot2lse, "Test", 1, emptyFilter, test_logger);
        
        inspector.insert(curLse);

        // All hits should be processed
        CHECK_FALSE(inspector.getSides().empty());
        CHECK(inspector.getSides().find(100)->second == 200);
    }

    SECTION("Active filter blocks missing hits") {
        std::map<int, std::set<int>> activeFilter;
        // Filter only knows query 999
        activeFilter[999] = {200};

        HitSetInspector inspector(expect, overlap, blastHits, prot2org, prot2len, allLse, prot2lse, "Test", 1, activeFilter, test_logger);
        
        // This should not throw, just skip the hit
        REQUIRE_NOTHROW(inspector.insert(curLse));

        // Result should be empty
        CHECK(inspector.getSides().empty());
    }

    SECTION("Filter is active, 100 is in it, but 200 is not allowed") {
        std::map<int, std::set<int>> activeFilter;
        // Only hits to subjects 300 and 400 are allowed for query 100
        activeFilter[100] = {300, 400};

        HitSetInspector inspector(expect, overlap, blastHits, prot2org, prot2len, allLse, prot2lse, "Test", 1, activeFilter, test_logger);
        
        inspector.insert(curLse);

        // Result empty as 200 is not in allowlist
        CHECK(inspector.getSides().empty());
    }

    SECTION("Filter is active, both proteins are allowed") {
        std::map<int, std::set<int>> activeFilter;
        activeFilter[100] = {200};

        HitSetInspector inspector(expect, overlap, blastHits, prot2org, prot2len, allLse, prot2lse, "Test", 1, activeFilter, test_logger);
        
        inspector.insert(curLse);

        // Hit has to be accepted
        CHECK(inspector.getSides().size() == 1);
        CHECK(inspector.getSides().find(100)->second == 200);
    }
}

TEST_CASE("HitSet::makeCOGs Cluster logic", "[cog][hitset]") {
    // Setup mock data
    DB db;
    
    db.prot2org[1] = "OrgA"; 
    db.prot2org[2] = "OrgA"; 
    db.prot2org[3] = "OrgB"; 
    db.prot2org[4] = "OrgB";
    db.prot2org[5] = "OrgC"; 

    db.allLse[1].lse = {1};
    db.allLse[2].lse = {2};
    db.allLse[3].lse = {3};
    db.allLse[4].lse = {4};
    db.allLse[5].lse = {5};

    // Dummy parameter
    double expect = 1e-10;
    double overlap = 0.5;
    std::string cName = "COG";
    int startNum = 1;
    
    std::map<int, std::set<int>> query2subject;
    auto test_logger = std::make_shared<spdlog::logger>("test", std::make_shared<spdlog::sinks::null_sink_st>());

    HitSetInspector inspector(
        expect, overlap, db.blastHits, db.prot2org, db.prot2len, 
        db.allLse, db.prot2lse, cName, startNum, query2subject, test_logger
    );

    SECTION("Error for edges between same organism") {
        auto& sides = const_cast<std::multimap<int, int>&>(inspector.getSides());
        // 1 and 2 are from OrgA
        sides.insert({1, 2});
        sides.insert({2, 1});

        CHECK_THROWS_AS(inspector.makeCOGs(), std::logic_error);
    }

    SECTION("Successful triangle detection across organisms") {
        auto& sides = const_cast<std::multimap<int, int>&>(inspector.getSides());
        
        sides.clear();
        sides.insert({1, 3}); sides.insert({3, 1}); // OrgA - OrgB
        sides.insert({1, 5}); sides.insert({5, 1}); // OrgA - OrgC
        sides.insert({3, 5}); sides.insert({5, 3}); // OrgB - OrgC

        COGResults results = inspector.makeCOGs();

        // 3 bidirectional edges
        CHECK(results.cog_edges.size() == 6);
        
        if (!results.cog_edges.empty()) {
            CHECK(results.cog_edges[0].cog_id.find("COG") != std::string::npos);
        }
    }

    SECTION("Detection of singletons") {
        // Register edge 1 to 3
        auto& sides = const_cast<std::multimap<int, int>&>(inspector.getSides());
        sides.clear();
        sides.insert({1, 3});
        sides.insert({3, 1});

        COGResults results = inspector.makeCOGs();

        bool found_protein_2_as_singleton = false;
    
    for (const auto& entry : results.cog_entries) {
        // We look for protein "2"
        if (entry.protein == "2") {
            // A singleton can be identified by the fact that cog_id is not set (std::nullopt)
            if (!entry.cog_id.has_value()) {
                found_protein_2_as_singleton = true;
            }
        }
    }
    
    CHECK(found_protein_2_as_singleton);
    }
}