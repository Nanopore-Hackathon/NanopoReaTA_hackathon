import groovy.json.JsonBuilder

process getParams {
    label "wf_common"
    cpus 1
    output:
        path "params.json"
    script:
        String paramsJSON = new JsonBuilder(params).toPrettyString() //Write down input parameters
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

