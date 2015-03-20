# Rakefile to collect all .tex files in a directory and run `pdflatex` and `bibtex` as needed to 
# produce PDF output.  If a .tex file is updated `pdflatex -draftmode` will be run to produce new 
# .aux and .log files.  These are used to determine whether `bibtex` needs to be run.  If so `bibtex` 
# will always need to be followed by `pdflatex -draftmode`.  With fully updated .aux and .bbl in 
# hand, a final `pdflatex` is run.  The only hole in the logic I've found is that, when making a 
# small revision, this will run `pdflatex -draftmode` then `pdflatex` when only `pdflatex` is 
# required.
#
# Run `rake` to compile PDFs and `rake clean` to remove the intermediary cruft

basedir = Dir.getwd 

TEX = FileList["**/*.tex"]
AUX = TEX.ext("aux")
BBL = TEX.ext("bbl")
BLG = TEX.ext("blg")
LOG = TEX.ext("log")
OUT = TEX.ext("out")
PDF = TEX.ext("pdf")

require 'rake/clean'
CLEAN.include(AUX)
CLEAN.include(BBL)
CLEAN.include(BLG)
CLEAN.include(LOG)
CLEAN.include(OUT)
CLOBBER.include(PDF)

desc "Full compile"
task :default => PDF 

desc "LaTeX aux"
rule ".aux" => ".tex" do |t|
	prefix = t.name.pathmap("%n")
	dir = t.name.pathmap("%d")
	Dir.chdir(dir)
	puts "pdflatex -draftmode #{prefix}"
	`pdflatex -draftmode #{prefix}`
	Dir.chdir(basedir)
end

desc "LaTeX log"
rule ".log" => ".tex" do |t|
	prefix = t.name.pathmap("%n")
	dir = t.name.pathmap("%d")
	Dir.chdir(dir)	
	puts "pdflatex -draftmode #{prefix}"
	`pdflatex -draftmode #{prefix}`
	Dir.chdir(basedir)	
end

desc "LaTeX compile"
# require log file and proceed if references are incomplete
rule ".pdf" => [".aux", ".bbl", ".log", ".tex"]  do |t|
	prefix = t.name.pathmap("%n")
	dir = t.name.pathmap("%d")
	Dir.chdir(dir)		
	puts "pdflatex #{prefix}"
	`pdflatex #{prefix}`
	Dir.chdir(basedir)		
end

desc "BibTex compile"
# look for .bib file in top-level directory
rule ".bbl" => [".aux", ".log", ".tex"]  do |t|
	prefix = t.name.pathmap("%n")
	dir = t.name.pathmap("%d")
	Dir.chdir(dir)			
	if cite?(prefix)
		puts "bibtex #{prefix}"
		`bibtex #{prefix}`
		puts "pdflatex -draftmode #{prefix}"
		`pdflatex -draftmode #{prefix}`
	end
	Dir.chdir(basedir)		
end

desc "Look at log file and check if references are complete"
def ref?(string)
	dirty = false
	file = File.open("#{string}.log", "r")
	m = file.read.match(/LaTeX Warning: There were undefined references|LaTeX Warning: Label(s) may have changed. Rerun to get|^LaTeX Warning: Reference/)
	file.close
	if m != nil
		puts m
		dirty = true
	end	
	return dirty
end

desc "Are citations up to date?"
def cite?(string)
	dirty = false
	if File.exists?("#{string}.bbl")
		if cite_missing(string) == true
			dirty = true
		else
			aux_list = cite_aux(string)
			bbl_list = cite_bbl(string)
			extra = (aux_list - bbl_list).length
			missing = (bbl_list - aux_list).length
			if extra > 0 || missing > 0
				dirty = true
			end
		end
	else
		dirty = true
	end
	return dirty
end

desc "Look at log file and check if citations are complete"
def cite_missing(string)
	dirty = false
	file = File.open("#{string}.log", "r")
	m = file.read.match(/^LaTeX Warning: Citation/)
	file.close
	if m != nil
		puts m
		dirty = true
	end	
	return dirty
end

desc "Find citations in .aux"
def cite_aux(string)
	list = Array.new
	file = File.open("#{string}.aux", "r")
	cites = file.read.scan(/\\citation{([^\}]+)}/)
	file.close
	if cites != nil
		cites.each {|m|
			list += m[0].split(",")
		}
	end	
	list.uniq!
	return list
end

desc "Find citations in .bbl"
def cite_bbl(string)
	list = Array.new
	file = File.open("#{string}.bbl", "r")
	cites = file.read.scan(/\\bibitem{([^\}]+)}/)
	file.close
	if cites != nil
		cites.each {|m|
			list += m[0].split(",")
		}
	end	
	list.uniq!
	return list
end