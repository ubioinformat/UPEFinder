This web is intended to help making queries to a Mongo database and displaying the results for everyone.
The Mongo database is provided in an external folder.

The app is created with Flask, and makes the queries with PyMongo

Dependencies:
    flask
    flask_caching
    bson
    flask_pymongo



Mongo database creation:
    1. start the mongo server: sudo service mongod start
    2. import the database from the json files: mongoimport --db=PRgene --collection=GeneCCLE --file=CCLE_Gene.json & mongoimport --db=PRgene --collection=GeneTcga --file=TCGA_Gene.json & mongoimport --db=PRgene --collection=GeneGtex --file=GTEX_Gene.json & mongoimport --db=PRgene --collection=Msig --file=Msig.json & mongoimport --db=PRgene --collection=GO --file=GO.json & mongoimport --db=PRgene --collection=Mala --file=Mala.json
    3. Enter mongo shell, select the database and create the indexes:
        mongo

        use PRgene

        db.GeneTcga.ensureIndex({"UPE":1 ,"GO_con.$id" : 1})
        db.GeneGtex.ensureIndex({"UPE":1 ,"GO_con.$id" : 1})
        db.GeneCCLE.ensureIndex({"UPE":1 ,"GO_con.$id" : 1})

        db.GeneTcga.ensureIndex({"UPE":1 ,"Mala_con.$id" : 1})
        db.GeneGtex.ensureIndex({"UPE":1 ,"Mala_con.$id" : 1})
        db.GeneCCLE.ensureIndex({"UPE":1 ,"Mala_con.$id" : 1})

        db.GeneTcga.ensureIndex({"UPE":1 ,"Msig_con.$id" : 1})
        db.GeneGtex.ensureIndex({"UPE":1 ,"Msig_con.$id" : 1})
        db.GeneCCLE.ensureIndex({"UPE":1 ,"Msig_con.$id" : 1})

Update routes:
In the run.py, two routes may be updated, the mongodb URI and the folder where the files to download will be created:
    app.config["MONGO_URI"] = 'mongodb://159.237.148.112:27017/PRgene'
    app.config["CLIENT_FILES"] = '/home/guille/data/FlaskTest2/tmp'


We planned to mount the server with nginX to serve the static files and with Gunicorn to serve the python script, with supervisor software to reset the server in case it crashes and to empty the tmp folder.The basic idea was to follow this tutorial: https://www.youtube.com/watch?v=goToXTC96Co