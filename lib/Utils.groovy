//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class Utils {

/*
NOTE: These can't be used in a config file to handle job submission in a better
way, but I'm just keeping them here for reference.

    // Extract numeric value and suffix from a string
    static def extract_number_and_suffix(value) {
        def match = (value =~ /([\d.]+)([a-zA-Z]+)/)
        if (match) {
            return [match[0][1] as double, match[0][2]]
        } else {
            throw new IllegalArgumentException("Invalid format: $value")
        }
    }

    // Adjust resources dynamically based on exit codes
    static def adjust_resources(base_value, resource_type, task) {
        if (task.attempt == 1) return base_value  // Keep defaults on first attempt

        def scale_factor = 1.0  // Initialize scale factor

        if ([71, 134, 137, 139, 140].contains(task.exitStatus) && resource_type == 'memory') {
            scale_factor = 1.5
        } else if ([143, 250].contains(task.exitStatus) && resource_type == 'time') {
            scale_factor = 1.5
        } else if ([104, 255].contains(task.exitStatus) && (resource_type == 'cpu' || resource_type == 'memory')) {
            scale_factor = 1.2
        }

        return base_value * scale_factor  // Return adjusted value
    }

    // Check if a queue has available CPU slots
    static def is_queue_available(queue_name, required_cpus) {
        def cmd = "qstat -g c | awk '\$1 == \"$queue_name\" {print \$3, \$4}'"
        def output = cmd.execute().text.trim()
        if (!output) return false
        def (used, total) = output.tokenize()*.toInteger()
        return (total - used) >= required_cpus
    }

    // Select the best queue based on CPU, memory, and time requirements
    static def select_queue(required_cpus, required_memory, required_time, params) {
        def (memory_value, memory_suffix) = extract_number_and_suffix(required_memory)
        def (time_value, time_suffix) = extract_number_and_suffix(required_time)
        for (queue in params.sge_queues.keySet()) {
            def resources = params.sge_queues[queue]
            def (queue_memory_value, _) = extract_number_and_suffix(resources.memory)
            def (queue_time_value, _unused_var) = extract_number_and_suffix(resources.time)
            if (resources.cpus >= required_cpus &&
                queue_memory_value >= memory_value &&
                queue_time_value >= time_value &&
                is_queue_available(queue, required_cpus)) {
                return queue
            }
        }
        return 'short.q'
    }
*/

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }
}
