import unittest
import json
from schema_salad.schema import load_schema
from schema_salad.validate import validate_ex
from schema_salad.sourceline import cmap

class TestValidate(unittest.TestCase):
    schema = cmap({"name": "_", "$graph":[{
        "name": "File",
        "type": "record",
        "fields": [{
            "name": "class",
            "type": {
                "type": "enum",
                "name": "File_class",
                "symbols": ["#_/File"]
            },
            "jsonldPredicate": {
                "_id": "@type",
                "_type": "@vocab"
            }
        }, {
            "name": "location",
            "type": "string",
            "jsonldPredicate": "_:location"
        }]
    }, {
        "name": "Directory",
        "type": "record",
        "fields": [{
            "name": "class",
            "type": {
                "type": "enum",
                "name": "Directory_class",
                "symbols": ["#_/Directory"]
            },
            "jsonldPredicate": {
                "_id": "@type",
                "_type": "@vocab"
            }
        }, {
            "name": "location",
            "type": "string",
            "jsonldPredicate": "_:location"
        }, {
            "name": "listing",
            "type": {
                "type": "array",
                "items": ["File", "Directory"]
            }
        }],
    }]})

    def test_validate_big(self):
        document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(self.schema)

        with open("biglisting.yml") as f:
            biglisting = json.load(f)

        self.assertEquals(True, validate_ex(avsc_names.get_name("Directory", ""), biglisting,
                                            strict=True, raise_ex=False))


    # def test_validate_small(self):
    #     document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(self.schema)

    #     with open("smalllisting.yml") as f:
    #         smalllisting = json.load(f)

    #     validate_ex(avsc_names.get_name("Directory", ""), smalllisting,
    #                 strict=True, raise_ex=True)
