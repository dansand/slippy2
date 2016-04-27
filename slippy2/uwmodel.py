#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

A class underworldModel that contains useful methods for pygplates - underworld integration

:license: GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""

import pygplates

import numpy as np
import pygplates

#class MyClass(object): = new-style class
#class MyClass: = OLD-STYLE CLASS

class underworldModel(object):
    """
    Class to create the necessary objects from a gPlates model (topologies, rotation_file)
    Main instance function is to restructure the model into age and feature type, as a dictionary.
    Other methods help with transformations.

    """
    def  __init__(self, topologies, rotation_file, timespan=(0,1,1), boundingBox=None):
        self.topologies = topologies
        self.rotation_file = rotation_file
        self.timespan = timespan
        self.boundingBox = boundingBox
        #Other class attrubutes
        self.coordsystem = "spherical"

    def get_times(self):
        """
        Returns a numpy array of the times (at which model features will be stored)
        from self.timespan = (start, stop, interval)
        """
        num = (self.timespan[0] - self.timespan[1])/ self.timespan[2]
        return np.linspace(self.timespan[1], self.timespan[0], (num+self.timespan[2]))

    def gplates_geoms(self, stage_pole_rotation=None, anchor_plate_id=0, fileformat = "gpml"):
        feature_dict = {}
        for age in self.get_times():
            ###########
            #Get topologies
            ###########
            ##########Not sure if the following should be instance attributes (self.resolved_topologies)
            time = age
            resolved_topologies = []
            resolved_topological_sections = []
            ############
            pygplates.resolve_topologies(self.topologies, self.rotation_file, resolved_topologies, time, resolved_topological_sections, anchor_plate_id=0)
            subfeats = []
            ############
            #Get unique boundary types
            ############
            feattypes = []
            for shared_boundary_section in resolved_topological_sections:
                for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            # Topological plate polygons and deforming networks have a boundary polygon with an area.
                    feat = shared_sub_segment.get_feature().get_feature_type()
                    featname = feat.get_name()
                    feattypes.append(featname)
            featurelist = list(np.unique((feattypes)))
            ############
            #Update dictionaries, and fill with dummy values
            Age = str(age)
            feature_type_dict = {name: 0 for name in featurelist}
            feature_dict.update({Age : feature_type_dict})
            for feat_name in featurelist:
                subfeats = []
                for shared_boundary_section in resolved_topological_sections:
                    for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                    # Topological plate polygons and deforming networks have a boundary polygon with an area.
                        if shared_sub_segment.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml(feat_name):
                            shared_sub_segment_feature = pygplates.Feature()
                            # Replace the cloned geometry property with the sub-segment geometry.
                            shared_sub_segment_feature.set_geometry(shared_sub_segment.get_geometry())
                            #Rotate all objects if required:
                            if stage_pole_rotation:
                                shared_sub_segment_feature.set_geometry(stage_pole_rotation*shared_sub_segment_feature.get_geometry())
                            subfeats.append(shared_sub_segment_feature)
                feature_type_dict.update({feat_name:subfeats})
        return feature_dict

    def uwmodel_transform(self,rotation=None, projection=None,  scale=None):
        """
        Might want to add something about spherical vs. cartesian
        if not projection:
        self.coordsystem = spherical
        """
        pass

    def bounding_box_transform(self, uwmodel_transform):
        pass
