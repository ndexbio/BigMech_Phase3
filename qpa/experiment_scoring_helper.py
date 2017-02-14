__author__ = 'aarongary'

from causal_paths.src.path_scoring import PathScoring
from os import path
import json

class ScoreExperiment:

    def __init__(self):
        self.place_holder = []

    def rank_score(self, paths):
        '''
        Ranks the provided paths according to the cross country algorithm and
        returns a dict in the format - {Target1: Rank, Target2: Rank, etc...}
        :param paths: dict of paths with Target as key
        :type paths: dictionary
        :return: dictionary of Target and path key value pairs
        :rtype: dict
        '''
        ps = PathScoring()

        score_these_path_tuples = []

        for i in sorted(paths):  # The path array i.e. [N1, E1, N2, E2, N3]
            replacement_path = []
            for j, path_item in enumerate(paths[i]):
                if j % 2 != 0:  # Odd elements are edges
                    replacement_path.append(self.convert_edge_dict_to_array(path_item))
                else:
                    replacement_path.append(path_item)

            path_tuple = ps.cx_edges_to_tuples(replacement_path, "A")

            #==========================================
            # Sum the individual edge scores to create
            # a total path score
            #==========================================
            sum = 0
            for tuple in path_tuple:
                sum += tuple[1]

            score_these_path_tuples.append((i, sum))

        scored_experiments = ps.calculate_average_position(score_these_path_tuples, [])

        #==========================================
        # Generates a dict with experiment Target
        # as key and rank float as value
        #==========================================
        edge_score_rank = {}

        for key in scored_experiments.keys():
            for e in scored_experiments[key]:
                edge_score_rank[e] = 1.0 / key

        for l in sorted(edge_score_rank):
            print "'%s': %f, " % (l, edge_score_rank[l])

        return edge_score_rank

    #==============================================
    # helper function to convert the raw edge dict
    # to an array which is the format used in
    # path scoring
    #==============================================
    def convert_edge_dict_to_array(self, edge):
        tmp_edge_list = []
        for e in edge.keys():

            tmp_edge_list.append(edge[e])

        return tmp_edge_list

'''
scoreExperiment = ScoreExperiment()

current_directory = path.abspath(path.dirname(__file__))

with open(path.join(current_directory, '../korkut_results.json'), 'r') as korkut_file:
    korkut_results = json.load(korkut_file)
    korkut_experiments = korkut_results.get("experiments")
    korkut_results_ranked = {}

    for exp_name in korkut_experiments.keys():
        try:
            target_paths = korkut_experiments.get(exp_name).get("target_paths")
            korkut_results_ranked[exp_name] = scoreExperiment.rank_score(target_paths)
        except KeyError:
            print "key error"

    print korkut_results_ranked

#scoreExperiment.rank_score(scoreExperiment.top_paths)

'''