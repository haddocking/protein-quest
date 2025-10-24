from protein_quest.alphafold.entry_summary import EntrySummary
from protein_quest.alphafold.fetch import converter


def test_loads_from_json_string():
    summary_json = '[{"toolUsed":"AlphaFold Monomer v2.0 pipeline","providerId":"GDM","entityType":"protein","isUniProt":true,"modelEntityId":"AF-P50613-F1","modelCreatedDate":"2025-08-01T00:00:00Z","sequenceVersionDate":"1996-10-01T00:00:00Z","globalMetricValue":82.0,"fractionPlddtVeryLow":0.162,"fractionPlddtLow":0.092,"fractionPlddtConfident":0.127,"fractionPlddtVeryHigh":0.618,"latestVersion":6,"allVersions":[1,2,3,4,5,6],"sequence":"MALDVKSRAKRYEKLDFLGEGQFATVYKARDKNTNQIVAIKKIKLGHRSEAKDGINRTALREIKLLQELSHPNIIGLLDAFGHKSNISLVFDFMETDLEVIIKDNSLVLTPSHIKAYMLMTLQGLEYLHQHWILHRDLKPNNLLLDENGVLKLADFGLAKSFGSPNRAYTHQVVTRWYRAPELLFGARMYGVGVDMWAVGCILAELLLRVPFLPGDSDLDQLTRIFETLGTPTEEQWPDMCSLPDYVTFKSFPGIPLHHIFSAAGDDLLDLIQGLFLFNPCARITATQALKMKYFSNRPGPTPGCQLPRPNCPVETLKEQSNPALAIKRKRTEALEQGGLPKKLIF","sequenceStart":1,"sequenceEnd":346,"sequenceChecksum":"0A94BFA7DD416CEB","isUniProtReviewed":true,"gene":"CDK7","uniprotAccession":"P50613","uniprotId":"CDK7_HUMAN","uniprotDescription":"Cyclin-dependent kinase 7","taxId":9606,"organismScientificName":"Homo sapiens","isUniProtReferenceProteome":true,"bcifUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.bcif","cifUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.cif","pdbUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-model_v6.pdb","paeImageUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-predicted_aligned_error_v6.png","msaUrl":"https://alphafold.ebi.ac.uk/files/msa/AF-P50613-F1-msa_v6.a3m","plddtDocUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-confidence_v6.json","paeDocUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-predicted_aligned_error_v6.json","amAnnotationsUrl":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-aa-substitutions.csv","amAnnotationsHg19Url":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-hg19.csv","amAnnotationsHg38Url":"https://alphafold.ebi.ac.uk/files/AF-P50613-F1-hg38.csv","entryId":"AF-P50613-F1","uniprotSequence":"MALDVKSRAKRYEKLDFLGEGQFATVYKARDKNTNQIVAIKKIKLGHRSEAKDGINRTALREIKLLQELSHPNIIGLLDAFGHKSNISLVFDFMETDLEVIIKDNSLVLTPSHIKAYMLMTLQGLEYLHQHWILHRDLKPNNLLLDENGVLKLADFGLAKSFGSPNRAYTHQVVTRWYRAPELLFGARMYGVGVDMWAVGCILAELLLRVPFLPGDSDLDQLTRIFETLGTPTEEQWPDMCSLPDYVTFKSFPGIPLHHIFSAAGDDLLDLIQGLFLFNPCARITATQALKMKYFSNRPGPTPGCQLPRPNCPVETLKEQSNPALAIKRKRTEALEQGGLPKKLIF","uniprotStart":1,"uniprotEnd":346,"isReferenceProteome":true,"isReviewed":true}]'

    results = converter.loads(summary_json, list[EntrySummary])

    assert isinstance(results, list)
    assert len(results) == 1
    assert all(isinstance(item, EntrySummary) for item in results)
    result = results[0]
    assert result.modelEntityId == "AF-P50613-F1"
    assert result.uniprotAccession == "P50613"
    assert result.fractionPlddtConfident == 0.127
