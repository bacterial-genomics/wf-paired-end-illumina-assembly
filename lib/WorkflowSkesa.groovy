//
// This file holds several functions specific to the workflow/skesa.nf in the wf-paired-end-illumina-assembly pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowSkesa {

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
