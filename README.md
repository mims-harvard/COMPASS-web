## Welcome to compass pages

## **compass 101**
#### This page provides details about the code, data, and pre-trained, fine-tuned model used for predicting responses to immunotherapy.


Date: 2024-03-17


## COMPASS 

[![Codebase](https://img.shields.io/badge/Codebase-Github-green)](https://github.com/mims-harvard/compass)
[![Document](https://img.shields.io/badge/Document-Github-yellow)](https://github.com/mims-harvard/compass-101/)
[![Slack](https://img.shields.io/badge/Project-Slack-orange)](https://zitniklab-harvard.slack.com/archives/C05S6LEQ3ED)


## Download all required files [here](https://drive.google.com/drive/folders/1ZsLvB9xUYHs4OrAXMXRYeKXNoaxMbMsG)

Please download all the files listed and save them to a folder on your local machine that you name, for example, the **./tmpignore** folder.

| File                                                                                                           | Description                                                                  |
| -------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| [compass](https://drive.google.com/open?id=1pcxzWXf-zq21TOjrCf1fHX6Dwgq8MThm&usp=drive_copy)                 | Conceptor main code, add to your python path to import                       |
| [pretrainer.pt](https://drive.google.com/open?id=1JnsKsacRiw0bRVpdP6DxbueFc0RmqYO9&usp=drive_copy)             | Pretrained model on TCGA data                                                |
| [finetuner_LFT_all.pt](https://drive.google.com/open?id=1dtI_y4R_FXA0SFKnUIp_--d4ufKJEgE8&usp=drive_copy)       | Fintuned model on ALL ITRP data (Linear probing)               |
| [finetuner_PFT_all.pt](https://drive.google.com/open?id=1fRWLtTJF49bSjO13kdC20Y0qtLYCvHkT&usp=drive_copy)       | Fintuned model on ALL ITRP data (Patial Finetuning)               |
| [ITRP_dataset.xlsx](https://drive.google.com/open?id=1oyIXRFrwRtC10gskbFmOffB54mW5O4MN&usp=drive_copy)         | The description of the ALL ITRP datasets                                     |
| [ITRP.PATIENT.TABLE](https://drive.google.com/open?id=1J5_b32-KhuSjxb8dOr9SVgG9ms3Bv21m&usp=drive_copy)        | The Patient clinical information of ITRP dataset                             |
| [ITRP.TPM.TABLE](https://drive.google.com/open?id=1VeKkzYjp51enwvFQw4OtxG08tfpYJYuT&usp=drive_copy)            | The Patient bulk-mRNA TPM value of ITRP dataset            

-----

## Install the dependency
In this example, we will demonstrate how to install the compass dependency and add it to our Python path for usage in this package.

* Step1: Go to the compass folder and install the dependency:
  ```bash
    cd compass
    pip install -r ./requirements.txt
  ```

* Step2: Now you can import conceptor by 

    ```python
    import sys
    sys.path.insert(0, './tmpignore/compass')
    from compass import FineTuner, loadcompass
    ```
