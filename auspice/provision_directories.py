import sys, os, shutil, re

# run this script to provision necessary directory structure and index.htmls
# this script deletes existing /h3n2/, /cdc/, etc... directories and makes new ones
# resulting directories are checked into git
# run from setup/ directory

def printdir():
    print(re.sub(r'^.+auspice/', 'auspice/', os.getcwd()))

indexfile = open("index.html", "w")
indexfile.write("---\n")
indexfile.write("title: nextflu\n")
indexfile.write("layout: redirect\n")
indexfile.write("rurl: /h3n2/ha/3y/\n")
indexfile.write("---\n")
indexfile.close()

# provision live site and gisaid site

virus_to_coloring = {
    "H3N2": "ep, ne, rb, lbi, dfreq, region, date, cHI",
    "H1N1pdm": "lbi, dfreq, region, date, cHI",
    "Vic": "lbi, dfreq, region, date, cHI",
    "Yam": "lbi, dfreq, region, date, cHI"
}

virus_to_freqdefault = {
    "H3N2": "1, 2, 3, 4, 5",
    "H1N1pdm": "6b.1, 6b.2, 164T",
    "Vic": "1A, 1B, 117V, DV",
    "Yam": "2, 3, 3a, 172Q"
}

viruses = ["H3N2", "H1N1pdm", "Vic", "Yam"]
segments = ["HA", "NA"]
resolutions = ["2y", "3y", "6y", "12y"]
for virus in viruses:

    vpath = virus.lower()
    if os.path.isdir(vpath):
        shutil.rmtree(vpath)
    os.makedirs(vpath)
    os.chdir(vpath)
    printdir()

    indexfile = open("index.html", "w")
    indexfile.write("---\n")
    indexfile.write("title: nextflu / %s \n" % virus)
    indexfile.write("layout: redirect\n")
    indexfile.write("rurl: /%s/ha/3y/\n" % vpath)
    indexfile.write("---\n")
    indexfile.close()

    for segment in segments:

        spath = segment.lower()
        if os.path.isdir(spath):
            shutil.rmtree(spath)
        os.makedirs(spath)
        os.chdir(spath)
        printdir()

        indexfile = open("index.html", "w")
        indexfile.write("---\n")
        indexfile.write("title: nextflu / %s / %s \n" % (virus, segment))
        indexfile.write("layout: redirect\n")
        indexfile.write("rurl: /%s/%s/3y/\n" % (vpath, spath))
        indexfile.write("---\n")
        indexfile.close()

        for resolution in resolutions:

            rpath = resolution.lower()
            if os.path.isdir(rpath):
                shutil.rmtree(rpath)
            os.makedirs(rpath)
            os.chdir(rpath)
            printdir()

            indexfile = open("index.html", "w")
            indexfile.write("---\n")
            indexfile.write("title: nextflu / %s / %s / %s\n" % (virus, segment, resolution))
            indexfile.write("layout: auspice\n")
            indexfile.write("virus: %s\n" % virus)
            indexfile.write("segment: %s\n" % segment)
            indexfile.write("resolution: %s\n" % resolution)
            indexfile.write("coloring: %s\n" % virus_to_coloring[virus])
            indexfile.write("gtplaceholder: HA1 positions...\n")
            indexfile.write("freqdefault: %s\n" % virus_to_freqdefault[virus])
            indexfile.write("site: public\n")
            indexfile.write("---\n")
            indexfile.write("\n")
            indexfile.write("<script>\n")
            indexfile.write("var file_prefix = \"flu_%s_%s_%s_\";\n" % (virus.lower(), segment.lower(), resolution.lower()))
            indexfile.write("var useTiters = false;\n")
            indexfile.write("{%% include %s_meta.js %%}\n" % virus)
            indexfile.write("{%% include %s_meta.js %%}\n" % resolution)
            indexfile.write("</script>\n")
            indexfile.close()

            os.chdir("..")

        os.chdir("..")

    for resolution in resolutions:

        rpath = resolution.lower()
        if os.path.isdir(rpath):
            shutil.rmtree(rpath)
        os.makedirs(rpath)
        os.chdir(rpath)
        printdir()

        indexfile = open("index.html", "w")
        indexfile.write("---\n")
        indexfile.write("title: nextflu / %s / %s \n" % (virus, resolution))
        indexfile.write("layout: redirect\n")
        indexfile.write("rurl: /%s/ha/%s/\n" % (vpath, rpath))
        indexfile.write("---\n")
        indexfile.close()

        os.chdir("..")

    os.chdir("..")

# provision GISAID site
if os.path.isdir("gisaid"):
    shutil.rmtree("gisaid")
os.makedirs("gisaid")
os.chdir("gisaid")
printdir()

indexfile = open("index.html", "w")
indexfile.write("---\n")
indexfile.write("title: nextflu\n")
indexfile.write("layout: redirect\n")
indexfile.write("rurl: /gisaid/h3n2/ha/3y/\n")
indexfile.write("---\n")
indexfile.close()

for virus in viruses:

    vpath = virus.lower()
    if os.path.isdir(vpath):
        shutil.rmtree(vpath)
    os.makedirs(vpath)
    os.chdir(vpath)
    printdir()

    indexfile = open("index.html", "w")
    indexfile.write("---\n")
    indexfile.write("title: nextflu / %s \n" % virus)
    indexfile.write("layout: redirect\n")
    indexfile.write("rurl: /gisaid/%s/ha/3y/\n" % vpath)
    indexfile.write("---\n")
    indexfile.close()

    for segment in segments:

        spath = segment.lower()
        if os.path.isdir(spath):
            shutil.rmtree(spath)
        os.makedirs(spath)
        os.chdir(spath)
        printdir()

        indexfile = open("index.html", "w")
        indexfile.write("---\n")
        indexfile.write("title: nextflu / %s / %s \n" % (virus, segment))
        indexfile.write("layout: redirect\n")
        indexfile.write("rurl: /gisaid/%s/%s/3y/\n" % (vpath, spath))
        indexfile.write("---\n")
        indexfile.close()

        for resolution in resolutions:

            rpath = resolution.lower()
            if os.path.isdir(rpath):
                shutil.rmtree(rpath)
            os.makedirs(rpath)
            os.chdir(rpath)
            printdir()

            indexfile = open("index.html", "w")
            indexfile.write("---\n")
            indexfile.write("title: nextflu / %s / %s / %s\n" % (virus, segment, resolution))
            indexfile.write("layout: auspice\n")
            indexfile.write("virus: %s\n" % virus)
            indexfile.write("segment: %s\n" % segment)
            indexfile.write("resolution: %s\n" % resolution)
            indexfile.write("coloring: %s\n" % virus_to_coloring[virus])
            indexfile.write("gtplaceholder: HA1 positions...\n")
            indexfile.write("freqdefault: %s\n" % virus_to_freqdefault[virus])
            indexfile.write("site: GISAID\n")
            indexfile.write("---\n")
            indexfile.write("\n")
            indexfile.write("<script>\n")
            indexfile.write("var file_prefix = \"flu_%s_%s_%s_\";\n" % (virus.lower(), segment.lower(), resolution.lower()))
            indexfile.write("var useTiters = false;\n")
            indexfile.write("{%% include %s_meta.js %%}\n" % virus)
            indexfile.write("{%% include %s_meta.js %%}\n" % resolution)
            indexfile.write("</script>\n")
            indexfile.close()

            os.chdir("..")

        os.chdir("..")

    for resolution in resolutions:

        rpath = resolution.lower()
        if os.path.isdir(rpath):
            shutil.rmtree(rpath)
        os.makedirs(rpath)
        os.chdir(rpath)
        printdir()

        indexfile = open("index.html", "w")
        indexfile.write("---\n")
        indexfile.write("title: nextflu / %s / %s \n" % (virus, resolution))
        indexfile.write("layout: redirect\n")
        indexfile.write("rurl: /gisaid/%s/ha/%s/\n" % (vpath, rpath))
        indexfile.write("---\n")
        indexfile.close()

        os.chdir("..")

    os.chdir("..")

os.chdir("..")

# provision WHO sites

virus_to_coloring = {
    "H3N2": "clade, age, age_score, glyc, ep, ne, rb, lbi, dfreq, region, date, cHI, HI_dist",
    "H1N1pdm": "age, age_score, glyc, lbi, dfreq, region, date, cHI, HI_dist",
    "Vic": "age, age_score, glyc, lbi, dfreq, region, date, cHI, HI_dist",
    "Yam": "age, age_score, glyc, lbi, dfreq, region, date, cHI, HI_dist"
}

virus_to_freqdefault = {
    "H3N2": "1, 2, 3, 4, 5",
    "H1N1pdm": "6b.1, 6b.2, 164T",
    "Vic": "1A, 1B, 117V, DV",
    "Yam": "2, 3, 3a, 172Q"
}

builds = ["CDC", "Crick", "NIID", "VIDRL", "WHO"]
viruses = ["H3N2", "H1N1pdm", "Vic", "Yam"]
segments = ["HA", "NA"]
resolutions = ["2y", "6y"]
passages = ["cell", "egg"]
assays = ["HI", "FRA"] # H3N2 only

for build in builds:

    bpath = build.lower()
    if os.path.isdir(bpath):
        shutil.rmtree(bpath)
    os.makedirs(bpath)
    os.chdir(bpath)
    printdir()

    indexfile = open("index.html", "w")
    indexfile.write("---\n")
    indexfile.write("title: nextflu / %s \n" % build)
    indexfile.write("layout: redirect\n")
    indexfile.write("rurl: /%s/h3n2/ha/2y/cell/hi/\n" % bpath)
    indexfile.write("---\n")
    indexfile.close()

    for virus in viruses:
        vpath = virus.lower()
        if os.path.isdir(vpath):
            shutil.rmtree(vpath)
        os.makedirs(vpath)
        os.chdir(vpath)
        printdir()

        indexfile = open("index.html", "w")
        indexfile.write("---\n")
        indexfile.write("title: nextflu / %s / %s \n" % (build, virus))
        indexfile.write("layout: redirect\n")
        indexfile.write("rurl: /%s/%s/ha/2y/cell/hi/\n" % (bpath, vpath))
        indexfile.write("---\n")
        indexfile.close()

        for segment in segments:

            spath = segment.lower()
            if os.path.isdir(spath):
                shutil.rmtree(spath)
            os.makedirs(spath)
            os.chdir(spath)
            printdir()

            indexfile = open("index.html", "w")
            indexfile.write("---\n")
            indexfile.write("title: nextflu / %s / %s / %s \n" % (build, virus, segment))
            indexfile.write("layout: redirect\n")
            indexfile.write("rurl: /%s/%s/%s/2y/cell/hi/\n" % (bpath, vpath, spath))
            indexfile.write("---\n")
            indexfile.close()

            for resolution in resolutions:

                rpath = resolution.lower()
                if os.path.isdir(rpath):
                    shutil.rmtree(rpath)
                os.makedirs(rpath)
                os.chdir(rpath)
                printdir()

                indexfile = open("index.html", "w")
                indexfile.write("---\n")
                indexfile.write("title: nextflu / %s / %s / %s / %s \n" % (build, virus, segment, resolution))
                indexfile.write("layout: redirect\n")
                indexfile.write("rurl: /%s/%s/%s/%s/cell/hi/\n" % (bpath, vpath, spath, rpath))
                indexfile.write("---\n")
                indexfile.close()

                for passage in passages:

                    ppath = passage.lower()
                    if os.path.isdir(ppath):
                        shutil.rmtree(ppath)
                    os.makedirs(ppath)
                    os.chdir(ppath)
                    printdir()

                    indexfile = open("index.html", "w")
                    indexfile.write("---\n")
                    indexfile.write("title: nextflu / %s / %s / %s / %s / %s \n" % (build, virus, segment, resolution, passage))
                    indexfile.write("layout: redirect\n")
                    indexfile.write("rurl: /%s/%s/%s/%s/%s/hi/\n" % (bpath, vpath, spath, rpath, ppath))
                    indexfile.write("---\n")
                    indexfile.close()

                    for assay in assays:

                        apath = assay.lower()
                        if os.path.isdir(apath):
                            shutil.rmtree(apath)
                        os.makedirs(apath)
                        os.chdir(apath)
                        printdir()

                        if virus != "H3N2" and assay == "FRA":

                            indexfile = open("index.html", "w")
                            indexfile.write("---\n")
                            indexfile.write("title: nextflu / %s / %s / %s / %s / %s / %s \n" % (build, virus, segment, resolution, passage, assay))
                            indexfile.write("layout: redirect\n")
                            indexfile.write("rurl: /%s/%s/%s/%s/%s/hi/\n" % (bpath, vpath, spath, rpath, ppath))
                            indexfile.write("---\n")
                            indexfile.close()

                        elif (virus == "Vic" or virus == "Yam") and build == "VIDRL":

                            indexfile = open("index.html", "w")
                            indexfile.write("---\n")
                            indexfile.write("title: nextflu / %s / %s / %s / %s / %s / %s \n" % (build, virus, segment, resolution, passage, assay))
                            indexfile.write("layout: redirect\n")
                            indexfile.write("rurl: /who/%s/%s/%s/%s/%s/\n" % (vpath, spath, rpath, ppath, apath))
                            indexfile.write("---\n")
                            indexfile.close()

                        else:

                            indexfile = open("index.html", "w")
                            indexfile.write("---\n")
                            indexfile.write("title: nextflu / %s / %s / %s / %s / %s / %s\n" % (build, virus, segment, resolution, passage, assay))
                            indexfile.write("layout: auspice\n")
                            indexfile.write("build: %s\n" % build)
                            indexfile.write("virus: %s\n" % virus)
                            indexfile.write("segment: %s\n" % segment)
                            indexfile.write("resolution: %s\n" % resolution)
                            indexfile.write("passage: %s\n" % passage)
                            indexfile.write("assay: %s\n" % assay)
                            indexfile.write("coloring: %s\n" % virus_to_coloring[virus])
                            indexfile.write("gtplaceholder: HA1 positions...\n")
                            indexfile.write("freqdefault: %s\n" % virus_to_freqdefault[virus])
                            indexfile.write("site: WHO\n")
                            indexfile.write("---\n")
                            indexfile.write("\n")
                            indexfile.write("<script>\n")
                            indexfile.write("var file_prefix = \"flu_%s_%s_%s_%s_%s_%s_\";\n" % (build.lower(), virus.lower(), segment.lower(), resolution.lower(), passage.lower(), assay.lower()))
                            indexfile.write("var useTiters = true;\n")
                            indexfile.write("{%% include %s_meta.js %%}\n" % virus)
                            indexfile.write("{%% include %s_meta.js %%}\n" % resolution)
                            indexfile.write("</script>\n")
                            indexfile.close()

                        os.chdir("..")

                    os.chdir("..")

                os.chdir("..")

            os.chdir("..")

        for resolution in resolutions:

            rpath = resolution.lower()
            if os.path.isdir(rpath):
                shutil.rmtree(rpath)
            os.makedirs(rpath)
            os.chdir(rpath)
            printdir()

            indexfile = open("index.html", "w")
            indexfile.write("---\n")
            indexfile.write("title: nextflu / %s / %s / %s \n" % (build, virus, resolution))
            indexfile.write("layout: redirect\n")
            indexfile.write("rurl: /%s/%s/ha/%s/cell/hi/\n" % (bpath, vpath, rpath))
            indexfile.write("---\n")
            indexfile.close()

            for passage in passages:

                ppath = passage.lower()
                if os.path.isdir(ppath):
                    shutil.rmtree(ppath)
                os.makedirs(ppath)
                os.chdir(ppath)
                printdir()

                indexfile = open("index.html", "w")
                indexfile.write("---\n")
                indexfile.write("title: nextflu / %s / %s / %s / %s \n" % (build, virus, resolution, passage))
                indexfile.write("layout: redirect\n")
                indexfile.write("rurl: /%s/%s/ha/%s/%s/hi/\n" % (bpath, vpath, rpath, ppath))
                indexfile.write("---\n")
                indexfile.close()

                for assay in assays:

                    if virus != "H3N2" and assay == "FRA":
                        continue

                    apath = assay.lower()
                    if os.path.isdir(apath):
                        shutil.rmtree(apath)
                    os.makedirs(apath)
                    os.chdir(apath)
                    printdir()

                    indexfile = open("index.html", "w")
                    indexfile.write("---\n")
                    indexfile.write("title: nextflu / %s / %s / %s / %s / %s \n" % (build, virus, resolution, passage, assay))
                    indexfile.write("layout: redirect\n")
                    indexfile.write("rurl: /%s/%s/ha/%s/%s/%s/\n" % (bpath, vpath, rpath, ppath, apath))
                    indexfile.write("---\n")
                    indexfile.close()

                    os.chdir("..")

                os.chdir("..")

            os.chdir("..")

        os.chdir("..")

    os.chdir("..")

# check
print("Directories provisioned")
