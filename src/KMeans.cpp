/** @file
* @Author: Andreas Dahl
* @Date:   2015-03-03 01:02:01
*/

#include "KMeans.h"

#include <algorithm>
#include <limits>

using std::vector;

static void find_centroid(vector<KMerString>& data, KMerString& best_kmer) {
    float best_distance = std::numeric_limits<float>::infinity();
    for (std::vector<KMerString>::iterator test = data.begin(); test != data.end(); ++test) {
        float current_dist = 0;
        for (std::vector<KMerString>::iterator test2 = data.begin(); test2 != data.end(); ++ test2) {
            current_dist += kMerDistanceHellinger(*test, *test);
        }
        if (current_dist < best_distance) {
            best_distance = current_dist;
            best_kmer = *test;
        }
    }
}

void kmeans(const vector<KMerString>& data, int k, vector<vector<KMerString>>& res) {
    // Initialize Centroids
    vector<KMerString> centroids;
    std::copy_n(data.begin(), k, std::back_inserter(centroids));
    while (true) {
        // Generate new clusters
        vector<vector<KMerString>> new_clusters(k);
        for (std::vector<KMerString>::const_iterator cur = data.begin(); cur != data.end(); ++cur) {
            // Find closest centroid
            float best_dist = std::numeric_limits<float>::infinity();
            int best_label = -1;
            for (unsigned int i = 0; i != centroids.size(); i++) {
                float tmp_dist = kMerDistanceHellinger(*cur, centroids[i]);
                if (tmp_dist < best_dist) {
                    best_dist = tmp_dist;
                    best_label = i;
                }
            }
            // add to centroid
            new_clusters[best_label].push_back(*cur);
        }
        // Check If iteration is same
        bool match = true;
        if (new_clusters.size() != res.size()) {
            match = false;
        }
        if (match) {
            for (unsigned int i = 0; i < new_clusters.size() && match; i++) {
                if (new_clusters[i].size() != res[i].size()) {
                    match = false;
                    break;
                }
                for (unsigned int j = 0; i < new_clusters[i].size(); i++) {
                    if (0 == kMerDistanceHellinger(new_clusters[i][j], res[i][j])) {
                        match = false;
                        break;
                    }
                }
            }
        }
        for (unsigned int i = 0; i < res.size(); ++i)
            res[i] = new_clusters[i];

        // React to match
        if (match) {
            return;
        }

        // Find new centroids
        for (unsigned int i = 0; i < res.size(); ++i) {
            find_centroid(res[i], centroids[i]);
        }
    }
}


