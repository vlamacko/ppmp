#!/usr/bin/env python3

from context import am

# Only run prediction.
all_df, single_df, triplet_df, all_test = am.analysis(whitney=False,
                                                      parallel=False,
                                                      dist=False,
                                                      violin=False,
                                                      order=False,
                                                      lvplot=False,
                                                      test_path='./data/test/',
                                                      sim=True)
