# from itertools import groupby
#
#
# class CrystalSitesList:
#     def __init__(self, crystal_sites=None):
#         if crystal_sites is None:
#             self._crystal_sites = []
#         else:
#             self._crystal_sites = crystal_sites
#
#         self._crystal_sites_by_index = {}
#
#         for site in self._crystal_sites:
#             self._crystal_sites_by_index[site.site_number].append(site)
#
#         self._mixed_site_indices = []
#         for index in self._crystal_sites_by_index.keys():
#             if len(self._crystal_sites_by_index[index]) > 1:
#                 self._mixed_site_indices.append(index)
#
#     def _reindex_sites(self):
#         """
#         Repopulates the crystal site and mixed site indices. Also renumbers the sites to consecutive integers.
#         :return:
#         """
#         self._crystal_sites_by_index = {}
#
#         for index, (_, sites) in enumerate(groupby(self._crystal_sites, lambda site: site.site_number)):
#             site_list = []
#             for site in sites:
#                 site.site_number = index
#                 site_list.append(site)
#             self._crystal_sites_by_index[index] = site_list
#
#         self._mixed_site_indices = []
#         for index in self._crystal_sites_by_index.keys():
#             if len(self._crystal_sites_by_index[index]) > 1:
#                 self._mixed_site_indices.append(index)
#
#     def add_site(self, new_site):
#         self._crystal_sites.append(new_site)
#         self._reindex_sites()
#
#     def remove_site(self, site_to_remove):
#         self._crystal_sites.remove(site_to_remove)
#         self._reindex_sites()
#
#     def to_string(self):
#         sorted_sites = sorted(self._crystal_sites, key=lambda site: site.site_number)
#         return "\n".join(map(lambda site: site.to_string(), sorted_sites))
#
#     def get_all_sites(self):
#         return self._crystal_sites
#
#     def get_sites_by_index(self, index):
#         return self._crystal_sites_by_index[index]
#
#     def get_mixed_site_indices(self):
#         return self._mixed_site_indices
#
#     def get_number_of_sites(self):
#         return len(self._crystal_sites_by_index)
