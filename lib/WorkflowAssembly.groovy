//
// This file holds several functions specific to the workflow/assembly.nf in the wf-paired-end-illumina-assembly pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowAssembly {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (!params.input) {
            log.error "Path to input samplesheet OR directory not specified."
            System.exit(1)
        }
    }
}
