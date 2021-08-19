#! /usr/bin/env python
# encoding: utf-8
# WARNING! Do not edit! https://waf.io/book/index.html#_obtaining_the_waf_file

import os,shutil
from waflib import Errors,Task,TaskGen,Utils,Node
@TaskGen.before_method('process_source')
@TaskGen.feature('buildcopy')
def make_buildcopy(self):
	def to_src_nodes(lst):
		if isinstance(lst,Node.Node):
			if not lst.is_src():
				raise Errors.WafError('buildcopy: node %s is not in src'%lst)
			if not os.path.isfile(lst.abspath()):
				raise Errors.WafError('buildcopy: Cannot copy directory %s (unsupported action)'%lst)
			return lst
		if isinstance(lst,str):
			lst=[x for x in Utils.split_path(lst)if x and x!='.']
		node=self.bld.path.get_src().search_node(lst)
		if node:
			if not os.path.isfile(node.abspath()):
				raise Errors.WafError('buildcopy: Cannot copy directory %s (unsupported action)'%node)
			return node
		node=self.bld.path.get_src().find_node(lst)
		if node:
			if not os.path.isfile(node.abspath()):
				raise Errors.WafError('buildcopy: Cannot copy directory %s (unsupported action)'%node)
			return node
		raise Errors.WafError('buildcopy: File not found in src: %s'%os.path.join(*lst))
	nodes=[to_src_nodes(n)for n in getattr(self,'buildcopy_source',getattr(self,'source',[]))]
	node_pairs=[(n,n.get_bld())for n in nodes]
	self.create_task('buildcopy',[n[0]for n in node_pairs],[n[1]for n in node_pairs],node_pairs=node_pairs)
class buildcopy(Task.Task):
	color='PINK'
	def keyword(self):
		return'Copying'
	def run(self):
		for f,t in self.node_pairs:
			t.parent.mkdir()
			shutil.copy2(f.abspath(),t.abspath())
