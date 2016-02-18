#!/usr/bin/python
# Severine Berard                                                Winter 2016
#
# Copyright 2015 Krister Swenson
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Stuff to deal with annotatons (under construction)
"""

from ncbi import fetchGBRecords, accFromDescription
from usefulfuncs import invdict
from Bio import Entrez


#        _________________________
#_______/        Functions          \_____________________________________________

def initializeIDTOREC(idTOstr):
  """
  Initialize the dictionnary IDTOREC

  @warn:   severine's code
  @param:  idTOstr, the dictionnary mapping genome id with their string description
  @return: the IDTOREC dictionnary
  """
  global IDTOREC
  idTOacc = {}
  Entrez.email = "severine.berard@umontpellier.fr" # getEmail() in sequencetool to be imported
  for genome in idTOstr:
      idTOacc[genome] = accFromDescription(idTOstr[genome])

  try:
    IDTOREC = fetchGBRecords(idTOacc.values(), invdict(idTOacc), './records')
  except IOError as e:
    sys.exit('Problem opening ./records directory.\n{}'.format(e))

def getAnnotations(gid,start,end):
  """
  Get the annotations recorded on genome gid between positions start and end

  @param gid:   genome id
  @param start: starting position of the subsequence of the genome to look for annotations
  @param end:   ending position  of the subsequence of the genome to look for annotations
  @param IDTOREC: the dictionnary (global variable ????)
  @return:      a list of annotations with start and end (annotation,start,end)
  """
  global IDTOREC
  information = []
  overlaprequired=0 #For the moment we take all overlapping annotations
  for f in IDTOREC[gid].features:
    s = int(f.location.start)
    e = int(f.location.end)
    if(('product' in f.qualifiers and (f.type == 'gene' or f.type == 'CDS'))
       and
       ((s >= start and e <= end) or # gene in interval
        (s < start and e > end)   or # interval in gene/CDS
        # 50% of annotation is in the interval
        ((e < end and (e-start)/float(e-s) > overlaprequired) or
         (s > start and (end-s)/float(e-s) > overlaprequired)))):
      annotation = f.qualifiers['product'][0]
      information.append((annotation,s,e))
  return information

