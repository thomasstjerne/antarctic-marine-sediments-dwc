const fs = require("fs");
const parse = require("csv-parse");
const transform = require("stream-transform");
const _ = require('lodash')
const Biom = require('biojs-io-biom').Biom;

const samples = require('./samples.json')
const mixs = require('./mixs.json')

const sampleByRun = _.keyBy(samples, 'run');
const mixsByID = _.keyBy(mixs, 'ID');
const PROJECT_ID = "PRJNA335729";
const VOTU_DB = "SILVA v132"
   const tax = {
       "sk": "kingdom",
       "p": "phylum",
       "c": "class",
       "o": "order",
       "f": "family",
       "g": "genus",
       "s": "species"
   }


  const getTaxon = (classification) => {
    let taxonomy = {
      "kingdom": "",
      "phylum": "",
      "class": "",
      "order": "",
      "family": "",
      "genus": ""
    };
    let scientificName = "";
    let rank = ""
    for(var i = 0; i < classification.length; i ++){
        const splitted = classification[i].split("__");
        if(splitted[1] && splitted[1] !== "" && tax[splitted[0]]){
            scientificName = splitted[1];
            rank = tax[splitted[0]];
            if(rank !== "species"){
              taxonomy[tax[splitted[0]]] = splitted[1];
            }
        }
    }
    return [taxonomy, scientificName, rank]

  }

const getMixs = (record, fastaMap) => {
    const seqRun = record.query.split('.')[0];
    const DNA_sequence = _.get(fastaMap, record.query, '');
    const sample = sampleByRun[seqRun];
    let silvaSplitted = record.SILVA.split(";");
    let i = silvaSplitted.length -1;
    while(i > -1 && silvaSplitted[i].endsWith("__")){
      silvaSplitted.pop();
      i--;
    }
    const occurrenceId = `${PROJECT_ID}:${sample.sample}:${silvaSplitted.join(";")}`;

    const {eventID} = sample;
    const sampleMixs = mixsByID[`${PROJECT_ID}__${eventID}`];
    const votu_db = VOTU_DB;
    const sample_name = sample.sample;
    const {project_name, env_broad_scale, env_local_scale, env_medium, experimental_factor, investigation_type, submitted_to_insdc, env_package, subspecf_gen_lin, pcr_primers, seq_meth, target_gene, lib_reads_seqd, url} = sampleMixs;
    const splitted_primers = pcr_primers.split(";");
    const pcr_primer_forward = splitted_primers[0].split(':')[1];
    const pcr_primer_reverse = splitted_primers[1].split(':')[1];
    const sop = "https://www.ebi.ac.uk/metagenomics/pipelines/4.1"
    const columns = [occurrenceId, submitted_to_insdc, `https://www.ebi.ac.uk/ena/browser/view/${project_name}`, `https://www.ebi.ac.uk/ena/browser/view/${sample_name}`, experimental_factor, env_broad_scale, env_local_scale, env_medium, env_package, subspecf_gen_lin, votu_db, target_gene, url, pcr_primers, pcr_primer_forward, pcr_primer_reverse, investigation_type, seq_meth, lib_reads_seqd, sop, DNA_sequence]
      return `${columns.join("\t")}\n`
};

const getEMOF = record => {
  const coreID =  record[0];
  const eventID = record[1];
  const sampleMixs = mixsByID[eventID];
  const {carb_nitro_ratio, chlorophyll, org_matter, phaeopigments, tot_diss_nitro, tot_org_carb} = sampleMixs;
  const columns = [
    {key: "carb_nitro_ratio", value: carb_nitro_ratio},
    {key: "chlorophyll", value: chlorophyll},
    {key: "org_matter", value: org_matter},
    {key: "phaeopigments", value: phaeopigments},
    {key: "tot_diss_nitro", value: tot_diss_nitro},
    {key: "tot_org_carb", value: tot_org_carb}
    ];
      return columns.map(c => `${coreID}\t${c.key}\t${c.value}`).join("\n")+"\n";
 
}

const writeOcc = (sampleID, classification, readCount, totalReads) => {
  const [taxonomy, scientificName, rank]  = getTaxon(classification);
  const sample = sampleByRun[sampleID];
  const {eventID} = sample;
    const sampleMixs = mixsByID[`${PROJECT_ID}__${eventID}`];
    const {lat_lon, collection_date, geo_loc_name} = sampleMixs;
  
    const splitLatLon = lat_lon.split(" ");
    const decimalLatitude = splitLatLon[0];
    const decimalLongitude = splitLatLon[1];
    const occurrenceId = `${PROJECT_ID}:${sample.sample}:${classification.join(';')}`;

        return `${occurrenceId}\t${PROJECT_ID}__${eventID}\t${collection_date}\t${geo_loc_name}\t${readCount}\t${totalReads}\tDNA sequence reads\tDNA sequence reads\t${decimalLatitude}\t${decimalLongitude}\t${scientificName}\t${Object.keys(taxonomy).map(rank => taxonomy[rank]).join("\t")}\t${rank}\thttps://www.ebi.ac.uk/ena/browser/view/${sample.sample}\tMATERIAL_SAMPLE`
 
}

const getFastaAsMap = (sample) => {
  const fasta = fs.readFileSync(`mgnify/${sample}_FASTQ_SSU.fasta`, 'utf8')
  const rows = fasta.split(">");
  let map = {};

  for(let i=0; i < rows.length -1; i++){
    if(rows[i]){
      const splitted = rows[i].split("\n");
      const ID = splitted[0];
      const sequence = splitted.slice(1).join("");
      map[ID] = sequence
    }
  }
  return map;
}
const processFile = async (sample) => {

 const biomString = fs.readFileSync(`mgnify/${sample}_FASTQ_SSU_OTU_TABLE_JSON.biom`, 'utf8')
  const biom = new Biom(JSON.parse(biomString));
  const readCounts = biom.getDataMatrix();
  const totalReads = readCounts.reduce((acc, curr) => acc+curr[0], 0)
  const taxa = biom.getMetadata({dimension: 'rows', attribute: 'taxonomy'});
  fs.writeFile(`processed/occurrence/${sample}_occ.txt`, taxa.map((t, index) => writeOcc(sample, t, readCounts[index][0], totalReads)).join("\n")+"\n", function (err) {
    if (err) return console.log(err);
    console.log(`occurrenceces written to ${sample}_occ.txt`);
  });
  const fastaMap = getFastaAsMap(sample);


  const parser = parse({
    delimiter: "\t",
    columns: true,
    ltrim: true,
    rtrim: true,
    quote: null
  });

  const occparser = parse({
    delimiter: "\t",
    columns: false,
    ltrim: true,
    rtrim: true,
    quote: null
  });


  const transformer = transform(function(record, callback){
      try {
        callback(null, getMixs(record, fastaMap))
      } catch(err){
          console.log(err)
          console.log(record)
      }
      
  }, {
    parallel: 5
  })

  const emofTransformer = transform(function(record, callback){
    try {
      callback(null, getEMOF(record))
    } catch(err){
        console.log(err)
        console.log(record)
    }
    
}, {
  parallel: 5
})
  var readStream = fs.createReadStream(__dirname+`/mgnify/${sample}_FASTQ_SSU_MAPSeq.mseq`);
  
  readStream.pipe(parser).pipe(transformer).pipe(fs.createWriteStream(__dirname+`/processed/sequence/${sample}_seq.txt`))
  
  var readStream2 = fs.createReadStream(__dirname+`/processed/occurrence/${sample}_occ.txt`);
  
  readStream2.pipe(occparser).pipe(emofTransformer).pipe(fs.createWriteStream(__dirname+`/processed/emof/${sample}_emof.txt`))

};

samples.forEach(s => processFile(s.run))
//processFile("SRR3997466");