#!/usr/bin/env python3

from context import am

# Run prediction and all the visualisation.
all_df, single_df, triplet_df, all_test = am.analysis(whitney=True,
                                                      parallel=True,
                                                      dist=True,
                                                      violin=True,
                                                      lvplot=True,
                                                      order=False,
                                                      test_path='./data/test/',
                                                      sim=True)
