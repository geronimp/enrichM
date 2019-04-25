#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


__author__ = "Ben Woodcroft, Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd3 near uq.edu.au"
__status__ = "Development"

###############################################################################
# Imports
import re
###############################################################################

KEGG = '^(K\d{5})$'
GH = '(GH\d+)'
PL = '(PL\d+)'
TIGRFAM = '^(TIGR\d{5})$'
PFAM = '(PF\d+)'
CE = '(CE\d+)'
EC = '\\d{1,2}(\\.(\\-|\\d{1,2})){3}'


class ModuleDescription:

    def __init__(self, module_description_string):
        self.module_description_string \
            = module_description_string
        self.parsed_module, self.is_single_step \
            = ModuleDescriptionParser().parse_module_string(module_description_string)

        if isinstance(self.parsed_module, ModuleDescriptionOrRelation):
            new_and = ModuleDescriptionAndRelation()
            new_and.relations = [self.parsed_module]
            self.parsed_module = new_and

    def kos(self):
        '''Return an iterable over the total list of KOs in the module'''
        r_kegg = re.compile(KEGG)
        r_gh = re.compile(GH)
        r_pl = re.compile(PL)
        r_tigrfam = re.compile(TIGRFAM)
        r_pfam = re.compile(PFAM)
        r_ce = re.compile(CE)
        r_ec = re.compile(EC)

        if any(r_kegg.findall(self.module_description_string)):
            return r_kegg.findall(self.module_description_string)
        elif any(r_gh.findall(self.module_description_string)):
            return r_gh.findall(self.module_description_string)
        elif any(r_pl.findall(self.module_description_string)):
            return r_pl.findall(self.module_description_string)
        elif any(r_tigrfam.findall(self.module_description_string)):
            return r_tigrfam.findall(self.module_description_string)
        elif any(r_pfam.findall(self.module_description_string)):
            return r_pfam.findall(self.module_description_string)
        elif any(r_ce.findall(self.module_description_string)):
            return r_ce.findall(self.module_description_string)
        elif any(r_ec.findall(self.module_description_string)):
            return r_ec.findall(self.module_description_string)

    def num_steps(self):
        if isinstance(self.parsed_module, ModuleDescriptionAndRelation):
            return len(self.parsed_module.relations)
        else:
            raise Exception("Cannot work with non-AND type modules")

    def num_covered_steps(self, ko_set):

        if isinstance(self.parsed_module, ModuleDescriptionAndRelation):
            step_cov = 0
            path_cov = 0
            reac_cov = 0
            ko_path = dict()

            for idx, m in enumerate(self.parsed_module.relations):
                step_passed, step_counts, reaction_counts, ko = m.satisfied_with(
                    ko_set, list())

                if step_passed:
                    step_cov += 1
                    path_cov += step_counts
                    ko_path[idx] = ko

                reac_cov += reaction_counts

            return step_cov, path_cov, reac_cov, ko_path

        else:
            raise Exception("Cannot work with non-AND type modules")


class GetLowestRelation:

    @staticmethod
    def get_lowest_relation(input_relation, to_fill):

        if type(input_relation) == ModuleDescriptionKoEntry:
            to_fill.append(input_relation.ko)
        else:

            for relation in input_relation.relations:
                if type(relation) == ModuleDescriptionKoEntry:
                    to_fill.append(relation.ko)
                else:
                    GetLowestRelation.get_lowest_relation(relation, to_fill)
        return to_fill


class ModuleDescriptionAndRelation:

    def satisfied_with(self, set_of_kos, kos):
        counts = 0
        step_passed = False
        reaction_counts = 0
        founds = list()

        for r in self.relations:
            found, count, reaction_count, ko = r.satisfied_with(
                set_of_kos, kos)

            if found:
                founds.append(1)
                counts += count

                for k in ko:

                    if k not in kos:
                        kos.append(k)

            reaction_counts += reaction_count

        step_passed = len(self.relations) == sum(founds)

        return step_passed, counts, reaction_counts, kos


class ModuleDescriptionPlusRelation(ModuleDescriptionAndRelation):
    pass


class ModuleDescriptionOrRelation:
    def satisfied_with(self, set_of_kos, kos):

        counts = 0
        step_passed = False
        reaction_counts = 0

        for r in self.relations:
            found, _, reaction_count, ko = r.satisfied_with(set_of_kos, kos)

            if found:
                step_passed = True
                for k in ko:
                    if k not in kos:
                        kos.append(k)
            else:
                kos_to_extract = GetLowestRelation.get_lowest_relation(
                    r, list())
                for ko in kos_to_extract:
                    if ko in kos:
                        kos.remove(ko)
            reaction_counts += reaction_count
        return step_passed, counts, reaction_counts, kos


class ModuleDescriptionKoEntry:
    def __init__(self, ko):
        self.ko = ko

    def satisfied_with(self, set_of_kos, kos):
        found = self.ko in set_of_kos
        count = (1 if found else 0)
        reaction_count = 1
        return found, count, reaction_count, [self.ko]


class ParserHelper:
    pass


class ModuleDescriptionParser:

    def correct_substrings(self, substring_list):
        fixed_substrings = list()
        for substring in substring_list:
            # Omit optional enzyme (e.g. M00372)
            # or undefined KO group (e.g. M00079)
            if substring.startswith('-'):
                continue
            # Remove redundant and definitions
            substring = substring.replace(', ', ',')
            fixed_substrings.append(substring)
        return fixed_substrings

    def parse_module_string(self, string):

        frags1 = self.split_on_space(string)
        frags1 = self.correct_substrings(frags1)
        is_single_step = len(frags1) == 1

        if is_single_step:
            if len(self.split_on_comma(frags1[0])) > 1:
                frags1 = self.split_on_comma(frags1[0])
                master_relation = ModuleDescriptionOrRelation()
            elif len(self.split_on_plus(frags1[0])) > 1:
                frags1 = self.split_on_plus(frags1[0])
                master_relation = ModuleDescriptionAndRelation()
            else:
                frags1 = self.split_on_comma(frags1[0])
                master_relation = ModuleDescriptionOrRelation()
        else:
            master_relation = ModuleDescriptionAndRelation()

        current = ParserHelper()
        current.top_relation = master_relation
        current.understuff = frags1
        stack = list([current])

        while len(stack) > 0:
            current = stack.pop()
            new_stuff = list()

            for e in current.understuff:

                if isinstance(e, str):

                    if (re.match(KEGG, e) or
                        re.match(GH, e) or
                        re.match(PL, e) or
                        re.match(CE, e) or
                        re.match(TIGRFAM, e) or
                        re.match(PFAM, e) or
                            re.match(EC, e)):

                        new_stuff.append(ModuleDescriptionKoEntry(e))

                    else:
                        frags = self.split_on_comma(e)

                        if len(frags) == 1:
                            topush = ParserHelper()
                            comma_splits = self.split_on_space(e)
                            m = None

                            if len(comma_splits) == 1:
                                plus_splits = self.split_on_plus(e)
                                minus_splits = self.split_on_minus(e)

                                if len(plus_splits) > 1:
                                    m = ModuleDescriptionPlusRelation()
                                    topush.understuff = plus_splits

                                elif len(minus_splits) > 1:
                                    m = ModuleDescriptionAndRelation()
                                    topush.understuff = minus_splits[:1]

                                else:
                                    raise Exception(
                                        "Parse exception on %s" % string)

                            else:
                                m = ModuleDescriptionAndRelation()
                                topush.understuff = comma_splits

                            topush.top_relation = m
                            stack.append(topush)
                            new_stuff.append(m)

                        else:
                            m = ModuleDescriptionOrRelation()
                            topush = ParserHelper()
                            topush.top_relation = m
                            topush.understuff = frags
                            stack.append(topush)
                            new_stuff.append(m)

                else:
                    new_stuff.append(e)

            current.top_relation.relations = new_stuff

        return master_relation, is_single_step

    def split_on(self, string, characters):
        bracket_counter = 0
        fragments = list()
        current = list()
        remove_end_brackets = True

        for i in range(len(string)):
            c = string[i]

            if c == '(':
                current += c
                bracket_counter += 1

            elif c == ')':
                current += c
                bracket_counter -= 1
                if bracket_counter == 0 and i < len(string)-1:
                    remove_end_brackets = False

            elif c in characters and bracket_counter == 0:
                fragments.append(''.join(current))
                current = ''

            else:
                current += c

        fragments.append(''.join(current))

        if remove_end_brackets and string[0] == '(':

            if string[-1] != ')':
                raise Exception("Parse error")
            return self.split_on(string[1:-1], characters)

        return fragments

    def split_on_space(self, string):
        return(self.split_on(string, ' '))

    def split_on_plus(self, string):
        return(self.split_on(string, '+'))

    def split_on_comma(self, string):
        return(self.split_on(string, ','))

    def split_on_minus(self, string):
        return(self.split_on(string, '-'))
