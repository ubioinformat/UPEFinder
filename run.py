import os, time
from flask import Flask, render_template, url_for, request, redirect, send_file, send_from_directory, session, safe_join
from flask_caching import Cache
from datetime import datetime
from flask import Blueprint
from flask_pymongo import PyMongo, ASCENDING, DESCENDING
from math import log2
from decimal import Decimal
import json
import uuid 

# cache = Cache(config={'CACHE_TYPE': 'simple'})

app = Flask(__name__)
app.config["MONGO_URI"] = 'mongodb://localhost:27017/PRgene'
app.config["CLIENT_FILES"] = os.path.join(os.getcwd(), 'tmp')
app.config["ANNOTATION_FILES"] = os.path.join(os.getcwd(), 'Annot')
app.secret_key = '123'
mongo = PyMongo(app)
# cache.init_app(app)


def get_bonds(mongoCursor, Ensg_set):
    results = []
    for query in mongoCursor:
        pe_conn = query['PE_con']
        pe_conn = set(map(lambda x: x['id'], pe_conn))
        intersect = Ensg_set.intersection(pe_conn)
        query['Number_bond'] = len(intersect)
        query['Names_bond']   = list(intersect)
        results.append(query)
    return results

def get_related_pval(relation_genes, target, database):
    if database == 'GO':
        pval_conn = 'GO_con'
    elif database == 'MSIG':
        pval_conn = 'Msig_con'
    elif database == 'MALA':
        pval_conn = 'Mala_con'
    for relation in range(0, len(relation_genes)):
        related_pval = list(filter(lambda conn: conn['id'] == target ,relation_genes[relation][pval_conn]))
        relation_genes[relation]['related_pval'] = '%.3E' % Decimal(float(related_pval[0]['pval']))
    return  sorted(relation_genes, key= lambda i: float(i['related_pval']), reverse = False)

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
    connections = {'GO' : "GO_con", 'MSIG' : "Msig_con", 'MALA' : "Mala_con"}
    collections = {'GO' : "GO", 'MSIG' : "Msig", 'MALA' : "Mala"}
    col = collections[database]
    connection =connections[database]
    htmlsimple = {'GO' : 'found_GO.html', 'MSIG' : 'found_Msig.html', 'MALA' : 'found_Mala.html'}
    htmldisamb = {'GO' : 'Disambiguation_GO.html', 'MSIG' : 'Disambiguation_MSIG.html', 'MALA' : 'Disambiguation_MALA.html'}
    try:
        print(field)
        print(target)
        stats = mongo.db[col].find({field : target}).next()
        print('Perfect match')
        if field != '_id':
            target_id = stats['_id']
            upesTCGA  = mongo.db.GeneTcga.find({"UPE" : True, connection : {'$elemMatch' : {'id':target_id}} })
            upesGTEX  = mongo.db.GeneGtex.find({"UPE" : True, connection : {'$elemMatch' : {'id':target_id}} })
            upesCCLE  = mongo.db.GeneCCLE.find({"UPE" : True, connection : {'$elemMatch' : {'id':target_id}} })
            pexgo_TCGA = set(map(lambda x: x['_id'], mongo.db.GeneTcga.find({"UPE" : False, connection : target_id }, {'_id':1})))
            pexgo_GTEX = set(map(lambda x: x['_id'], mongo.db.GeneGtex.find({"UPE" : False, connection : target_id }, {'_id':1})))
            pexgo_CCLE = set(map(lambda x: x['_id'], mongo.db.GeneCCLE.find({"UPE" : False, connection : target_id }, {'_id':1})))
        else:
            upesTCGA = mongo.db.GeneTcga.find({"UPE" : True, connection : {'$elemMatch' : {'id':target}} })
            upesGTEX = mongo.db.GeneGtex.find({"UPE" : True, connection : {'$elemMatch' : {'id':target}} })
            upesCCLE = mongo.db.GeneCCLE.find({"UPE" : True, connection : {'$elemMatch' : {'id':target}} })
            pexgo_TCGA = set(map(lambda x: x['_id'], mongo.db.GeneTcga.find({"UPE" : False, connection : {'$elemMatch' : {'id':target}} }, {'_id':1})))
            pexgo_GTEX = set(map(lambda x: x['_id'], mongo.db.GeneGtex.find({"UPE" : False, connection : {'$elemMatch' : {'id':target}} }, {'_id':1})))
            pexgo_CCLE = set(map(lambda x: x['_id'], mongo.db.GeneCCLE.find({"UPE" : False, connection : {'$elemMatch' : {'id':target}} }, {'_id':1})))
        genesXGO_TCGA = get_bonds(upesTCGA, pexgo_TCGA)
        genesXGO_TCGA = get_related_pval(genesXGO_TCGA, target, database)
        genesXGO_GTEX = get_bonds(upesGTEX, pexgo_GTEX)
        genesXGO_GTEX = get_related_pval(genesXGO_GTEX, target, database)
        genesXGO_CCLE = get_bonds(upesCCLE, pexgo_CCLE)
        genesXGO_CCLE = get_related_pval(genesXGO_CCLE, target, database)
        return htmlsimple[database], stats, genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target
    except StopIteration:
        print('pattern match')
        stats = mongo.db[col].find({field : {'$regex' : target, '$options': 'i'}})
        return htmldisamb[database], stats, None, None, None, target
    return

def search_for_download_with_bonds(database, upe_id, target, field):
    databases = {'TCGA' : 'GeneTcga', 'GTEX' : 'GeneGtex', 'CCLE' : 'GeneCCLE'}
    fields = {'GO' : 'GO_con.id', 'MSIG' : 'Msig_con.id', 'MALA' : 'Mala_con.id'}
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
# @cache.cached(timeout=50)
def find_Ensg(ID, target):
    template, stats, gene_tcga, gene_gtex, gene_ccle, target = simple_search('_id', ID)
    if (gene_tcga.count() !=0) or (gene_tcga.count() !=0) or (gene_tcga.count() !=0):
        # Perfect match
        print('Perfect match')
        return render_template(template, stats=stats, genesTCGA=gene_tcga, genesGTEX=gene_gtex, genesCCLE=gene_ccle, target=target)
    else:
        return '<h1> SORRY!!</h1>'


@app.route('/find_GO/<target>/<ID>', methods = ['GET', 'POST'])
# @cache.cached(timeout=50)
def find_Go(ID, target):
    template, go_stats , genesXGO_TCGA, genesXGO_GTEX, genesXGO_CCLE, target = complex_search('_id', ID, 'GO')
    return render_template(template, go_stats=go_stats, genesXGO_TCGA=genesXGO_TCGA, genesXGO_GTEX=genesXGO_GTEX, genesXGO_CCLE=genesXGO_CCLE, target=target)

@app.route('/find_MSIG/<target>/<ID>', methods = ['GET', 'POST'])
def find_Msig(ID, target):
    template, Msig_stats , genesXMSIG_TCGA, genesXMSIG_GTEX, genesXMSIG_CCLE, target = complex_search('_id', ID, 'MSIG')
    return render_template(template, msig_stats=Msig_stats, genesXMSIG_TCGA=genesXMSIG_TCGA, genesXMSIG_GTEX=genesXMSIG_GTEX, genesXMSIG_CCLE=genesXMSIG_CCLE, target=target)


@app.route('/find_MALA/<target>/<ID>', methods = ['GET', 'POST'])
# @cache.cached(timeout=50)
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
        HEADER = ['\t'.join(['Dataset', 'neXtProt_id', 'gene_id', 'gene_name', 'uPE1', 'Connections_with_PE1/uPE1', 'GO_connections', 'MSigdb_connections', 'Disease_connections', 'N_PE_bonding_toquery'])+'\n']
        genesX_TCGA = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_TCGA)
        genesX_GTEX = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_GTEX)
        genesX_CCLE = map(lambda x: SEP.join([x['dataset'], x['nexprot_id'], x['_id'], x['gene_name'],str( x['UPE']), str(x['N_PE_con']), str(x['N_GO_con']), str(x['N_MSIG_con']), str(x['N_Mala_con']), str(x['Number_bond'])],) , genesX_CCLE)
    else:
        _, _, genesX_TCGA, genesX_GTEX, genesX_CCLE, _ = simple_search('_id', target)
        HEADER = ['\t'.join(['Dataset', 'neXtProt_id', 'gene_id', 'gene_name', 'uPE1', 'Connections_with_PE1/uPE1', 'GO_connections', 'MSigdb_connections', 'Disease_connections'])+'\n']
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
    entry = search_for_download_with_bonds(database=database, upe_id=ID, target=target, field=origin)[0]
    gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
    gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
    header = ['Dataset', 'NexprotID', 'EnsemblID', 'GeneName', 'uPE1', 'PE_level', 'Chr', 'Bonding', 'Coexpressed_EnsemblID', 'Correlation',\
            'Coexpressed_NextprotID', 'Coexpressed_GeneName', 'Coexpressed_gene_type', \
            'Coexpressed_Chr', 'Coexpressed_start', 'Coexpressed_end',	'Coexpressed_ProteinEvidence',\
            'Coexpressed_PE','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Coexpr_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
        PE_connections = entry['PE_con']
        bonding_genes = entry['Names_bond']
        for PE in PE_connections:
            PE_name = PE['id']
            if PE_name in bonding_genes:
                bonding = '1'
            else:
                bonding = '0'
            line = anchor +SEP+ bonding +SEP+ PE_name +SEP+ str(PE['corr']) +SEP+  SEP.join(gene_annot[PE_name][1:3])+SEP+  SEP.join(gene_annot[PE_name][4:8])+'\n'
            outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

@app.route('/download_coexpressed/<target>/<database>/<ID>', methods = ['GET', 'POST'])
def download_coexpressed(target, database, ID):
    entry = search_for_download(database, ID).next()
    gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
    gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
    header = ['Dataset', 'NexprotID', 'EnsemblID', 'GeneName', 'uPE1', 'PE_level', 'Chr', 'Coexpressed_EnsemblID', \
        'Coexpressed_NextprotID', 'Coexpressed_GeneName', 'Correlation', 'Coexpressed_gene_type', \
        'Coexpressed_Chr', 'Coexpressed_start', 'Coexpressed_end','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Coexpr_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
        PE_connections = entry['PE_con']
        for PE in PE_connections:
            line = anchor +SEP+ PE['id'] +SEP+  SEP.join(gene_annot[PE['id']][1:3]) +SEP+ str(PE['corr']) +SEP+  SEP.join(gene_annot[PE['id']][4:8])+'\n'
            outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

@app.route('/download_Annotations/<target>/<database>/<ID>', methods = ['GET', 'POST'])
def download_Annotations(target, database, ID):
    print(target, database, ID)
    entry = search_for_download(database, ID).next()
    annot = open(os.path.join(app.config['ANNOTATION_FILES'],'GOannotation_Ago20.txt'), 'r').readlines()
    del annot[0]
    annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), annot))
    header = ['Dataset', 'neXtProt_ID', 'EnsembleID', 'gene_name', 'uPE1', 'PE_level', 'Chr', 'DB', 'Ontology', 'P.value','Description','\n'] if entry['UPE'] else ['Dataset', 'neXtProt_ID', 'EnsembleID', 'gene_name', 'uPE1', 'PE_level', 'Chr', 'DB', 'Ontology', 'Description','\n']
    SEP = '\t'
    file_name = target.replace(' ','_').title() + 'Annot_' +time.strftime("%Y%m%d-%H%M") + '.tsv'
    with open(os.path.join('tmp', file_name), 'w') as outfl:
        outfl.write(SEP.join(header))
        anchor = SEP.join([entry['dataset'], entry['nexprot_id'], entry['_id'], entry['gene_name'], str(entry['UPE']), entry['PE_level'], str(entry['chr'])])
        for go in entry['GO_con']:
            go_id   = go['id']
            if entry['UPE']:
                go_pval = go['pval']
                line = anchor +SEP+ 'Gene Ontology' +SEP+ go_id +SEP+ go_pval +SEP+ annot[go_id][2] + '\n'
            else:
                line = anchor +SEP+ 'Gene Ontology' +SEP+ go_id +SEP+ annot[go_id][2] + '\n'
            outfl.write(line)
        for msig in entry['Msig_con']:
            msig_id   = msig['id']
            if entry['UPE']:
                msig_pval = msig['pval']
                line = anchor +SEP+ 'MSigDB' +SEP+ msig_id +SEP+ msig_pval +'\n'
            else:
                line = anchor +SEP+ 'MSigDB' +SEP+ msig_id +'\n'
            outfl.write(line)
        for mala in entry['Mala_con']:
            mala_id   = mala['id']
            if entry['UPE']:
                mala_pval = mala['pval']
                line = anchor +SEP+ 'Disease' +SEP+ mala_id +SEP+ mala_pval + '\n'
            else:
                line = anchor +SEP+ 'Disease' +SEP+ mala_id + '\n'
            outfl.write(line)
    safe_path = safe_join(app.config["CLIENT_FILES"], file_name)
    return send_file( safe_path, as_attachment=True)

def get_PRscore(collection, database, ID):
    collections = {'GO' : "GO", 'MSIG' : "Msig", 'MALA' : "Mala", 'TCGA': "GeneTcga", 'GTEX': "GeneGtex", 'CCLE': "GeneCCLE"}
    col = collections[collection]
    fields = {'TCGA':"gene_name", 'CCLE':"gene_name", 'GTEX':"gene_name", 'GO':"Term", 'MSIG':'_id', 'MALA': '_id'}
    fld = fields[collection]
    scores = {'TCGA': 'PR_score', 'GTEX': 'PR_score', 'CCLE': 'PR_score',}
    if collection not in ('TCGA', 'GTEX', 'CCLE'):
        pr_score = mongo.db[col].find({fld:ID}, {'_id' : 0 , 'PR_score' + '_' + database : 1}).next()
        return log2(float(pr_score['PR_score' + '_' + database]))
    else:
        pr_score = mongo.db[col].find({fld:ID}, {'_id' : 0 , 'PR_score' : 1}).next()
        if not pr_score['PR_score']:
            return float(-15)
        return log2(float(pr_score['PR_score']))

def generate_template_cytoscape_breadFirst(ID, database):
    _ , _, gene_tcga, gene_gtex, gene_ccle, _ = simple_search('_id', ID)
    # only for TCGA at the moment
    if database == 'TCGA':
        genes = gene_tcga.next()
    elif database == 'GTEX':
        genes = gene_gtex.next()
    elif database == 'CCLE':
        genes = gene_ccle.next()
    root_type = 'UPE' if genes['UPE'] else 'Gene'
    node_type = 'UPE' if not genes['UPE'] else 'Gene'
    # Load annotation files 
    gene_annot = open(os.path.join(app.config['ANNOTATION_FILES'], 'nextprot_all_genecode+unknown.txt'), 'r').readlines()
    gene_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), gene_annot))
    go_annot = open(os.path.join(app.config['ANNOTATION_FILES'],'GOannotation_Ago20.txt'), 'r').readlines()
    del go_annot[0]
    go_annot = dict(map(lambda x: (x.split('\t')[0], x.split('\t')), go_annot))
    # PE nodes
    pe_nodes = list(map(lambda x: x['id'], genes['PE_con']))
    pe_nodes_ann = list(map(lambda x: gene_annot[x][2], pe_nodes))
    go_nodes = list(map(lambda x: x['id'], genes['GO_con']))
    go_nodes_ann =  list(map(lambda x: go_annot[x][2], go_nodes))
    msig_nodes = list(map(lambda x: x['id'], genes['Msig_con']))
    mala_nodes = list(map(lambda x: x['id'], genes['Mala_con']))
    nodes_list = list(map(lambda x: {'data':{'id':x, 'type': node_type, 'PR_score' : get_PRscore(database, database, x),'corr' : next(filter( lambda conn: conn['gene_name'] == x ,genes['PE_con']))['corr']}}, pe_nodes_ann))
    nodes_list = nodes_list + list(map(lambda x: {'data':{'id':x, 'type': 'GO', 'PR_score' : get_PRscore('GO', database, x)}}, go_nodes_ann)) 
    nodes_list = nodes_list + list(map(lambda x: {'data':{'id':x, 'type': 'Msig', 'PR_score' : get_PRscore('MSIG', database, x)}}, msig_nodes))
    nodes_list = nodes_list + list(map(lambda x: {'data':{'id':x, 'type': 'Mala', 'PR_score' : get_PRscore('MALA', database, x)}}, mala_nodes))
    # Get the pr ranges before adding the central node
    pr_ranges = list(map(lambda x: x['data']['PR_score'], nodes_list))
    corr_ranges = list(map(lambda x: x['corr'],genes['PE_con']))
    # add the UPE
    nodes_list.append({'data':{'id':genes['gene_name'], 'type': root_type, 'PR_score' : 0, 'corr' : 1, 'center':'True'}})
    # edges to the UPE
    all_nodes = pe_nodes_ann #+go_nodes_ann+msig_nodes+mala_nodes
    # all_nodes = pe_nodes_ann + go_nodes_ann + msig_nodes + mala_nodes
    edges = list(map(lambda x: {'data':{'source': genes['gene_name'], 'target':x , 'type' : 'Gene'}} , all_nodes))
    # edges = edges + list(map(lambda x: {'data':{'source': x, 'target':genes['_id'] , 'type' : 'Gene'}} , all_nodes))
    # get all the edges
    collections = {'GO' : "GO", 'MSIG' : "Msig", 'MALA' : "Mala", 'TCGA': "GeneTcga", 'GTEX': "GeneGtex", 'CCLE': "GeneCCLE"}
    col = collections[database]
    go_edges = []
    msig_edges = []
    mala_edges = []
    for pe in pe_nodes:
        conn = mongo.db[col].find({ '_id':pe}).next()
        go_conn = list(map(lambda x: x['id'], conn['GO_con']))
        msig_conn = list(map(lambda x: x['id'], conn['Msig_con']))
        mala_conn = list(map(lambda x: x['id'], conn['Mala_con']))
        intersection = list(set(go_conn) & set(go_nodes))
        go_edges.extend(list(map(lambda x: {'data':{'source': go_annot[x][2], 'target': gene_annot[pe][2], 'type': 'GO'}} , intersection)))
        intersection = list(set(msig_conn) & set(msig_nodes))
        msig_edges.extend(list(map(lambda x: {'data':{'source': x, 'target': gene_annot[pe][2], 'type': 'Msig'}} , intersection)))
        intersection = list(set(mala_conn) & set(mala_nodes))
        msig_edges.extend(list(map(lambda x: {'data':{'source': x, 'target': gene_annot[pe][2], 'type': 'Mala'}} , intersection)))
    edges = edges + go_edges + msig_edges + mala_edges
    # unique identifier
    unique_ID = uuid.uuid1().int
    # json to file
    json_name = f'json{unique_ID}.json'
    json_data = json.dumps(nodes_list + edges, indent=4, sort_keys=True)
    # print(safe_join(app.config["CLIENT_FILES"], json_name))
    # with open(safe_join(app.config["CLIENT_FILES"], json_name), 'w') as outfl:
    #     outfl.write(json_data) 
    # with open(safe_join(app.config["CLIENT_FILES"], f'json{unique_ID}.js'), 'w') as outfl:
    #     outfl.write('var upe_data = {\n')
    #     outfl.write('"nodes":')
    #     tmp1 = json.dumps(nodes_list, indent=4, sort_keys=True)
    #     outfl.write(tmp1[:-1])
    #     outfl.write(',], edges:')
    #     tmp1 = json.dumps(edges, indent=4, sort_keys=True)
    #     outfl.write(tmp1)
    #     outfl.write('};')
    infl =  open('./templates/cytoscape_template_custom_BreadFirst.html', 'r')
    print(unique_ID)
    safe_name = f'cytoscape_custom_{unique_ID}.html'
    print(safe_name)
    safe_path = safe_join(app.config["CLIENT_FILES"], safe_name)
    print(safe_path)
    new_template = open(safe_path, 'w')
    for line in infl:
        if 'JSONDATA' in line:
            newline = line.replace('JSONDATA', json_data)
            new_template.write(newline)
        elif ('MIN_pr' in line) and ('MAX_pr' in line):
            newline = line.replace('MIN_pr', str(round(min(pr_ranges), 3))).replace('MAX_pr', str(round(max(pr_ranges),3)))
            new_template.write(newline)
        elif ('MIN_corr' in line) and ('MAX_corr' in line):
            newline = line.replace('MIN_corr', str(round(min(corr_ranges), 3))).replace('MAX_corr', str(round(max(corr_ranges),3)))
            new_template.write(newline)
        else: 
            new_template.write(line)
    new_template.close()
    return safe_path 

@app.route('/get_graph/<database>/<ID>', methods = ['GET', 'POST'])
def get_graph(ID, database):
    print('graph')
    print(ID)
    safe_path = generate_template_cytoscape_breadFirst(ID, database)
    return  send_file(safe_path, as_attachment=False)

if __name__ == '__main__':
    app.run(debug = True, port=5000,  host='0.0.0.0')




