"""Train models"""

import click as c


@c.command
@c.option(
    "--model",

    "model_choice",
    help="Which model to use",
    default="model1",
    show_default=True,
    type=c.STRING, 
)
def train(
    model_choice,
): 
    if model_choice == "model1":
        raise NotImplementedError
    elif model_choice == "model2":

        model 



if __name__ == "__main__":

    train()