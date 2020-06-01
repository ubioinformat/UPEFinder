import os, time
from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, session, safe_join
from flask_caching import Cache
from datetime import datetime
from flask import Blueprint
from flask_pymongo import PyMongo, ASCENDING, DESCENDING
from itertools import chain
from bson import json_util


cache = Cache(config={'CACHE_TYPE': 'simple'})

app = Flask(__name__)
app.config["MONGO_URI"] = 'mongodb://localhost:27017/PRgene'
app.config["CLIENT_FILES"] = os.path.join(os.getcwd(), 'tmp')
app.config["ANNOTATION_FILES"] = os.path.join(os.getcwd(), 'Annot')
app.secret_key = '123'
mongo = PyMongo(app)
cache.init_app(app)


def get_bonds(mongoCursor, Ensg_set):
    results = []
    for query in mongoCursor:
        pe_conn = query['PE_con']
        pe_conn = set(map(lambda x: x.id, pe_conn))
        intersect = Ensg_set.intersection(pe_conn)
        query['Number_bond'] = len(intersect)
        query['Names_bond']   = list(intersect)
        results.append(query)
    return results


def simple_search(field, target):
    print('simple_search')
    gene_tcga   = mongo.db.GeneTcga.find({field : target}).sort([(field, ASCENDING)])
    N_gene_tcga = mongo.db.GeneTcga.count_documents({field : target})
    gene_gtex   = mongo.db.GeneGtex.find({field : target}).sort([(field, ASCENDING)])
    N_gene_gtex = mongo.db.GeneGtex.count_documents({field : target})
    gene_ccle   = mongo.db.GeneCCLE.find({field : target}).sort([(field, ASCENDING)])
    N_gene_ccle = mongo.db.GeneCCLE.count_documents({field : target})
    if (N_gene_tcga !=0) or (N_gene_gtex !=0) or (N_gene_ccle !=0):
        # Perfect match
        print('Perfect match')
        try:
            stats = mongo.db.GeneTcga.find({field : target}).next()
        except StopIteration:
            stats = mongo.db.GeneGtex.find({field : target}).next()
        gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
        gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
        stats = dict(zip(['ENSG', 'NextprotID', 'GeneName', 'gene_id', 'gene_type', 'Chr', 'start', 'end', 'ProteinEvidence', 'PE', 'Missing', 'strand', 'unknown_function'], gene_annot[stats['_id']]))
        return 'found_ENSG.html', stats, gene_tcga, gene_gtex, gene_ccle, target
    else:
        #chek for pattern match
        print('pattern match')
        gene_tcga = mongo.db.GeneTcga.find({field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        gene_gtex = mongo.db.GeneGtex.find({field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        gene_ccle = mongo.db.GeneCCLE.find({field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        return 'Disambiguation_ENSG.html', None, gene_tcga, gene_gtex, gene_ccle, target

def simple_search_upe(field, target):
    gene_tcga = mongo.db.GeneTcga.find({'UPE':True, field : target}).sort([(field, ASCENDING)])
    N_gene_tcga = mongo.db.GeneTcga.count_documents({'UPE':True, field : target})
    gene_gtex = mongo.db.GeneGtex.find({'UPE':True, field : target}).sort([(field, ASCENDING)])
    N_gene_gtex = mongo.db.GeneGtex.count_documents({'UPE':True, field : target})
    gene_ccle = mongo.db.GeneCCLE.find({'UPE':True, field : target}).sort([(field, ASCENDING)])
    N_gene_ccle = mongo.db.GeneCCLE.count_documents({'UPE':True, field : target})
    if (N_gene_tcga !=0) or (N_gene_gtex !=0) or (N_gene_ccle !=0):
        # Perfect match
        print('Perfect match')
        return 'found_ENSG.html', gene_tcga, gene_gtex, gene_ccle, target
    else:
        #chek for pattern match
        print('pattern match')
        gene_tcga = mongo.db.GeneTcga.find({'UPE':True, field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        gene_gtex = mongo.db.GeneGtex.find({'UPE':True, field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        gene_ccle = mongo.db.GeneCCLE.find({'UPE':True, field:{'$regex': target, '$options': 'i'}}).sort([(field, ASCENDING)])
        return 'Disambiguation_ENSG.html', gene_tcga,gene_gtex, gene_ccle, target


def complex_search(field, target, database):
    connections = {'GO' : "GO_con.$id", 'MSIG' : "Msig_con.$id", 'MALA' : "Mala_con.$id"}
    collections = {'GO' : "GO", 'MSIG' : "Msig", 'MALA' : "Mala"}
    col = collections[database]
    connection =connections[database]
    htmlsimple = {'GO' : 'found_GO.html', 'MSIG' : 'found_Msig.html', 'MALA' : 'found_Mala.html'}
    htmldisamb = {'GO' : 'Disambiguation_GO.html', 'MSIG' : 'Disambiguation_MSIG.html', 'MALA' : 'Disambiguation_MALA.html'}
    try:
        print(field)
        print(mongo.db[col])
        print(target)
        stats = mongo.db[col].find({field : target}).next()
        print('Perfect match')
        if field != '_id':
            target_id = stats['_id']
            upesTCGA  = mongo.db.GeneTcga.find({"UPE" : True, connection : target_id })
            upesGTEX  = mongo.db.GeneGtex.find({"UPE" : True, connection : target_id })
            upesCCLE  = mongo.db.GeneCCLE.find({"UPE" : True, connection : target_id })
            pexgo_TCGA = set(map(lambda x: x['_id'], mongo.db.GeneTcga.find({"UPE" : False, connection : target_id }, {'_id':1})))
            pexgo_GTEX = set(map(lambda x: x['_id'], mongo.db.GeneGtex.find({"UPE" : False, connection : target_id }, {'_id':1})))
            pexgo_CCLE = set(map(lambda x: x['_id'], mongo.db.GeneCCLE.find({"UPE" : False, connection : target_id }, {'_id':1})))
        else:
            upesTCGA = mongo.db.GeneTcga.find({"UPE" : True, connection : target })
            upesGTEX = mongo.db.GeneGtex.find({"UPE" : True, connection : target })
            upesCCLE = mongo.db.GeneCCLE.find({"UPE" : True, connection : target })
            pexgo_TCGA = set(map(lambda x: x['_id'], mongo.db.GeneTcga.find({"UPE" : False, connection : target }, {'_id':1})))
            pexgo_GTEX = set(map(lambda x: x['_id'], mongo.db.GeneGtex.find({"UPE" : False, connection : target }, {'_id':1})))
            pexgo_CCLE = set(map(lambda x: x['_id'], mongo.db.GeneCCLE.find({"UPE" : False, connection : target }, {'_id':1})))
        genesXGO_TCGA = get_bonds(upesTCGA, pexgo_TCGA)
        genesXGO_GTEX = get_bonds(upesGTEX, pexgo_GTEX)
        genesXGO_CCLE = get_bonds(upesCCLE, pexgo_CCLE)
        return htmlsimple[database], stats, genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target
    except StopIteration:
        print('pattern match')
        stats = mongo.db[col].find({field : {'$regex' : target, '$options': 'i'}})
        return htmldisamb[database], stats, None, None, None, target
    return

def search_for_download_with_bonds(database, upe_id, target, field):
    databases = {'TCGA' : 'GeneTcga', 'GTEX' : 'GeneGtex', 'CCLE' : 'GeneCCLE'}
    fields = {'GO' : 'GO_con.$id', 'MSIG' : 'Msig_con.$id', 'MALA' : 'Mala_con.$id'}
    upes = mongo.db[databases[database]].find({"_id" : upe_id})
    pes = set(map(lambda x: x['_id'], mongo.db[databases[database]].find({"UPE" : False, fields[field] : target }, {'_id':1})))
    gene_with_bonds= get_bonds(upes, pes)
    return gene_with_bonds

def search_for_download(database, upe_id):
    databases = {'TCGA' : 'GeneTcga', 'GTEX' : 'GeneGtex', 'CCLE' : 'GeneCCLE'}
    genes = mongo.db[databases[database]].find({"_id" : upe_id})
    return genes

@app.route('/', methods = ['POST', 'GET'])
def index():
    if request.method == 'POST':
        target = request.form['content'].strip()
        print(target)
        # if no ID to search print some demo
        if not target:
            gene_tcga = mongo.db.GeneTcga.find().limit(10)
            gene_gtex = mongo.db.GeneGtex.find().limit(10)
            gene_ccle = mongo.db.GeneCCLE.find().limit(10)
            return render_template('found_ENSG.html', genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle)
        ###############
        ##CLIENT_INPUTS
        ###############
        print(request.form.get('drop-down'))
        if request.form.get('drop-down') == 'neXtProt':
            print('neXtProt selected')
            template, stats, gene_tcga, gene_gtex, gene_ccle, target = simple_search('nexprot_id', target)
            return render_template(template, stats=stats, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'Ensemble':
            print('Ensemble selected')
            template, stats, gene_tcga, gene_gtex, gene_ccle, target= simple_search('_id', target)
            return render_template(template, stats=stats, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'Gene_Name':
            print('Gene_name selected')
            template, stats, gene_tcga, gene_gtex, gene_ccle, target = simple_search('gene_name', target)
            return render_template(template, stats=stats, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'GO_Term':
            print('GO_Term selected')
            template, go_stats , genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target = complex_search('Term', target, 'GO')
            return render_template(template, go_stats=go_stats, genesXGO_TCGA=genesXGO_TCGA, genesXGO_GTEX=genesXGO_GTEX, genesXGO_CCLE=genesXGO_CCLE, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'GO_Id':
            print('GO_ID selected')
            template, go_stats , genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target = complex_search('_id', target, 'GO')
            print(target)
            return render_template(template, go_stats=go_stats, genesXGO_TCGA=genesXGO_TCGA, genesXGO_GTEX=genesXGO_GTEX, genesXGO_CCLE=genesXGO_CCLE, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'MSIG':
            print('MSIG selected')
            template, Msig_stats , genesXMSIG_TCGA, genesXMSIG_GTEX, genesXMSIG_CCLE, target = complex_search('Description', target, 'MSIG')
            return render_template(template, msig_stats=Msig_stats, genesXMSIG_TCGA=genesXMSIG_TCGA, genesXMSIG_GTEX=genesXMSIG_GTEX, genesXMSIG_CCLE=genesXMSIG_CCLE, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'Disease':
            print('Disease selected')
            template, Mala_stats , genesXMALA_TCGA, genesXMALA_GTEX, genesXMALA_CCLE, target = complex_search('_id', target, 'MALA')
            return render_template(template, mala_stats=Mala_stats, genesXMALA_TCGA=genesXMALA_TCGA, genesXMALA_GTEX=genesXMALA_GTEX, genesXMALA_CCLE=genesXMALA_CCLE, target=target, choice=request.form.get('drop-down'))

        if request.form.get('drop-down') == 'Chromosome':
            print('Chromosome selected')
            template, gene_tcga, gene_gtex, gene_ccle, target = simple_search_upe('chr', target)
            return render_template(template, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target, choice=request.form.get('drop-down'))

    else:
        # Start Page
        return render_template('Welcome.html')



@app.route('/find_ENSG/<target>/<ID>', methods = ['GET', 'POST'])
@cache.cached(timeout=50)
def find_Ensg(ID, target):
    template, stats, gene_tcga, gene_gtex, gene_ccle, target = simple_search('_id', ID)
    if (gene_tcga.count() !=0) or (gene_tcga.count() !=0) or (gene_tcga.count() !=0):
        # Perfect match
        print('Perfect match')
        return render_template(template, stats=stats, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target)
    else:
        return '<h1> SORRY!!</h1>'


@app.route('/find_GO/<target>/<ID>', methods = ['GET', 'POST'])
@cache.cached(timeout=50)
def find_Go(ID, target):
    template, go_stats , genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target = complex_search('_id', ID, 'GO')
    return render_template(template, go_stats=go_stats, genesXGO_TCGA=genesXGO_TCGA, genesXGO_GTEX=genesXGO_GTEX, genesXGO_CCLE=genesXGO_CCLE, target=target)

@app.route('/find_MSIG/<target>/<ID>', methods = ['GET', 'POST'])
def find_Msig(ID, target):
    template, Msig_stats , genesXMSIG_TCGA, genesXMSIG_GTEX, genesXMSIG_CCLE, target = complex_search('_id', ID, 'MSIG')
    return render_template(template, msig_stats=Msig_stats, genesXMSIG_TCGA=genesXMSIG_TCGA, genesXMSIG_GTEX=genesXMSIG_GTEX, genesXMSIG_CCLE=genesXMSIG_CCLE, target=target)


@app.route('/find_MALA/<target>/<ID>', methods = ['GET', 'POST'])
@cache.cached(timeout=50)
def find_Mala(ID, target):
    template, Mala_stats , genesXMALA_TCGA, genesXMALA_GTEX, genesXMALA_CCLE, target = complex_search('_id', ID, 'MALA')
    return render_template(template, mala_stats=Mala_stats, genesXMALA_TCGA=genesXMALA_TCGA, genesXMALA_GTEX=genesXMALA_GTEX, genesXMALA_CCLE=genesXMALA_CCLE, target=target)


@app.route('/download_all/<target>/<field>', methods = ['GET', 'POST'])
def download_all(target, field):
    print(target, field)
    SEP = '\t'

    if field == 'GO':
        _, _, genesX_TCGA, genesX_GTEX, genesX_CCLE, _ = complex_search('_id', target, 'GO')
    if field == 'MSIG':
        _, _, genesX_TCGA, genesX_GTEX, genesX_CCLE, _ = complex_search('_id', target, 'MSIG')
    if field == 'MALA':
        _, _, genesX_TCGA, genesX_GTEX, genesX_CCLE, _ = complex_search('_id', target, 'MALA')
    if  field != 'ENSG':
        HEADER = ['\t'.join(['dataset', 'nexprot_id', '_id', 'gene_name', 'uPE1', 'Connections_with_PE1/uPE1', 'GO_connections', 'MSigdb_connections', 'Disease_connections', 'N_PE_bonding_toquery'])+'\n']
        genesX_TCGA = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_TCGA)
        genesX_GTEX = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_GTEX)
        genesX_CCLE = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_CCLE)
    else:
        _, _, genesX_TCGA, genesX_GTEX, genesX_CCLE, _ = simple_search('_id', target)
        HEADER = ['\t'.join(['dataset', 'nexprot_id', '_id', 'gene_name', 'uPE1', 'Connections_with_PE1/uPE1', 'GO_connections', 'MSigdb_connections', 'Disease_connections'])+'\n']
        genesX_TCGA = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con'])]) , genesX_TCGA)
        genesX_GTEX = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con'])]) , genesX_GTEX)
        genesX_CCLE = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con'])]) , genesX_CCLE)
    target = target.replace(' ', '_').title()
    file_name = target + '_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp',file_name), 'w') as outfl:
        outfl.writelines(HEADER)
        outfl.writelines('\n'.join(genesX_TCGA)+'\n')
        outfl.writelines('\n'.join(genesX_GTEX)+'\n')
        outfl.writelines('\n'.join(genesX_CCLE))
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    try:
        return send_file( safe_path, as_attachment=True)
    except:
        return '<h1>SORRY!!</h1>'


@app.route('/download_coexpressed_PEgenes/<origin>/<target>/<database>/<ID>', methods = ['GET', 'POST'])
def download_coexpressed_PEgenes(origin, target, database, ID):
    print(target, database, ID)
    data = search_for_download_with_bonds(database=database, upe_id=ID, target=target, field=origin)
    gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
    gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
    header = ['Dataset', 'NexprotID', 'EnsemblID', 'GeneName', 'uPE1', 'PE_level', 'Chr', 'Bonding', 'Coexpressed_EnsmbleID', \
            'Coexpressed_NextprotID', 'Coexpressed_GeneName', 'Coexpressed_gene_type', \
            'Coexpressed_Chr', 'Coexpressed_start', 'Coexpressed_end',	'Coexpressed_ProteinEvidence',\
            'Coexpressed_PE','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Coexpr_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        for entry in data:
            anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
            PE_connections = list(map(lambda x: x.id, entry['PE_con']))
            bonding_genes = entry['Names_bond']
            for PE_name in PE_connections:
                if PE_name in bonding_genes:
                    bonding = '1'
                else:
                    bonding = '0'
                line = anchor +SEP+ bonding +SEP+ PE_name +SEP+  SEP.join(gene_annot[PE_name][1:3])+SEP+  SEP.join(gene_annot[PE_name][4:8])+'\n'
                outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

@app.route('/download_coexpressed/<target>/<database>/<ID>', methods = ['GET', 'POST'])
def download_coexpressed(target, database, ID):
    data = search_for_download(database, ID)
    gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
    gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
    header = ['Dataset', 'NexprotID', 'EnsemblID', 'GeneName', 'uPE1', 'PE_level', 'Chr', 'Coexpressed_EnsmbleID', \
        'Coexpressed_NextprotID', 'Coexpressed_GeneName', 'Coexpressed_gene_type', \
        'Coexpressed_Chr', 'Coexpressed_start', 'Coexpressed_end','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Coexpr_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        for entry in data:
            anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
            PE_connections = list(map(lambda x: x.id, entry['PE_con']))
            for PE_name in PE_connections:
                line = anchor +SEP+ PE_name +SEP+  SEP.join(gene_annot[PE_name][1:3])+SEP+  SEP.join(gene_annot[PE_name][4:8])+'\n'
                outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

@app.route('/download_Annotations/<target>/<database>/<ID>', methods = ['GET', 'POST'])
def download_Annotations(target, database, ID):
    print(target, database, ID)
    data = search_for_download(database, ID)
    annot = open(os.path.join(app.config['ANNOTATION_FILES'],'GOannotation.txt'), 'r').readlines()
    del annot[0]
    annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), annot))
    header = ['Dataset', 'NexprotID', 'EnsmbleID', 'GeneName', 'uPE1', 'PE_level', 'Chr', 'DB', 'Ontology', 'Description','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Annot_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        for entry in data:
            anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
            GO_connections = list(map(lambda x: x.id, entry['GO_con']))
            for go in GO_connections:
                line = anchor +SEP+ 'Gene Ontology' +SEP+ go +SEP+ annot[go][2] + '\n'
                outfl.write(line)
            Msig_connections = list(map(lambda x: x.id, entry['Msig_con']))
            for msig in Msig_connections:
                line = anchor +SEP+ 'MSigDB' +SEP+ msig+'\n'
                outfl.write(line)
            Mala_connections = list(map(lambda x: x.id, entry['Mala_con']))
            for mala in Mala_connections:
                line = anchor +SEP+ 'Disease' +SEP+ mala+'\n'
                outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

if __name__ == '__main__':
    app.run(debug = True, port=5000,  host='0.0.0.0')
