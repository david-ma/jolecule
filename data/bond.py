import math
import json

import pdbstruct
import vector3d
from spacehash import SpaceHash


def get_bonds(atoms):
  separation_sq = lambda a, b: vector3d.pos_distance_sq(a.pos, b.pos)
  atom_to_vertex = lambda a: [a.pos.x, a.pos.y, a.pos.z]
  CHONPS = ['C', 'H', 'O', 'N', 'P', 'S']
  small_cutoff = 1.2**2
  medium_cutoff = 1.9**2
  large_cutoff = 2.4**2
  bonds = []
  space_hash = SpaceHash(map(atom_to_vertex, atoms))
  for i_atom, j_atom in space_hash.close_pairs():
    a = atoms[i_atom]
    b = atoms[j_atom]
    if (a.element == "H") or (b.element == "H"):
      cutoff = small_cutoff
    elif (a.element in CHONPS) and (b.element in CHONPS):
      cutoff = medium_cutoff
    else: 
      # probably a metal
      cutoff = large_cutoff
    if separation_sq(a, b) < cutoff:
      bonds.append([i_atom, j_atom])
  extent = math.sqrt(sum(map(lambda s:s**2, space_hash.spans)))
  return bonds, extent


def extract_bonds_max_length(pdb, js_fname):
  polymer = pdbstruct.Polymer(pdb)
  atoms = polymer.atoms()
  bonds, extent = get_bonds(atoms)  
  f = open(js_fname, 'w')
  lines = [l for l in open(pdb).readlines()
           if l.startswith("ATOM") or l.startswith("HETATM")]
  f.write('var lines = %s;' % json.dumps(lines))
  f.write('\n\n')
  f.write('var bond_pairs = %s;' % json.dumps(bonds))
  f.write('\n\n')
  f.write('var max_length = %s;' % json.dumps(extent))
  f.write('\n\n')
  f.write('var filename = "";')
  f.close()


def convert_js_to_json(js_fname, json_fname):
  js_text = open(js_fname).read()
  i_lines = js_text.index('var lines =')
  i_bonds = js_text.index('var bond_pairs =')
  i_max_length = js_text.index('var max_length =')
  i_filename = js_text.index('var filename = ')
  lines = js_text[i_lines:i_bonds]
  bonds = js_text[i_bonds:i_max_length]
  max_length = js_text[i_max_length:i_filename]
  filename = js_text[i_filename:-1]

  def strip_replace(s, substring):
    result = s.replace(substring, '').strip()
    if result[-1] == ';':
      result = result[:-1]
    return result
  lines = strip_replace(lines, 'var lines =')
  bonds = strip_replace(bonds, 'var bond_pairs =')
  max_length = strip_replace(max_length, 'var max_length =')
  filename = strip_replace(filename, 'var filename =')
  protein_data = {
    'pdb_atom_lines': eval(lines),
    'bond_pairs': eval(bonds),
    'max_length': eval(max_length),
    'filename': eval(filename),
  }
  json_str = 'var protein_data = ' + json.dumps(protein_data)
  open(json_fname, 'w').write(json_str)



extract_bonds_max_length('1mbo.pdb', '1mbo.data.js')

convert_js_to_json('1mbo.data.js', '1mbo.data.json')

  