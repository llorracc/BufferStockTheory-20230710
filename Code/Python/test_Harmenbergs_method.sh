#!/bin/bash

# Install requirements for Harmenberg test, if not already installed
if [ ! -e requirements_Harmenberg.out ]; then
    export PYTHONPATH=$(python -c "import site, os; print(os.path.join(site.USER_SITE,'lib','python','site-packages'))"):$PYTHONPATH:./HarKmenberg/
    [[ ! -e HarKmenberg ]] && mkdir HarKmenberg
    pip install --target HarKmenberg -r requirements_Harmenberg.txt | tee requirements_Harmenberg.out
fi

# Descend to the point in the HARK toolkit where the tests are executed
pushd .

# Execute only the test of Harmenberg's method
python -m pytest --log-cli-level=DEBUG src/econ-ark/HARK/ConsumptionSaving/tests/test_IndShockConsumerType.py -k test_Harmenbergs_method

# Return to original point
popd 

