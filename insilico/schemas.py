from schema import (Optional, Or, Schema)


# Specification of a task to apply a filter
schema_apply_filter = Schema({
    "apply_filter": {
        "id": int,
        Optional("depends_on"): int,
        "filters": {
            "functional_groups": {
                "smiles": Schema([str])
            }
        }
    }
}
)

# Specification of a task to compute a property
schema_compute_property = Schema({
    "compute_property": {
        "id": int,
        Optional("depends_on"): int,
        "property": Schema([str])
    }
}
)

# Specification of a task to search for a property
schema_search_property = Schema({
    "search_property": {
        "id": int,
        Optional("depends_on"): int,
        "property": Schema([str])
    }
}
)


schema_workflow = Schema({
    "name": str,
    "file_molecules": str,
    "steps": Schema([Or(schema_apply_filter, schema_compute_property, schema_search_property)])
})
