"""YAML parser for configuration object"""
from functools import partial
from pathlib import Path
from typing import Callable, Union

import yaml
from pydantic import BaseModel

from CASM.protein.config import (
    ClusterConfig,
)

def config_constructor(
    loader: yaml.FullLoader, node: yaml.nodes.MappingNode
) -> BaseModel:
    """
    Construct a BaseModel config.
    
    :param loader: Given YAML loader
    :param type: yaml.FullLoader
    :param node: A mapping node
    :param type: yaml.nodes.MappingNode
    """
    arg_map = loader.construct_mapping(node, deep=True) if node.value else {}
    return eval(node.tag[1:])(**arg_map)

def function_constructor(
    loader: yaml.FullLoader,
    tag_suffix: str,
    node: Union[yaml.nodes.MappingNode, yaml.nodes.ScalarNode],
) -> Callable:
    """
    Construct a Callable.  

    :param loader: Given YAML loader
    :param type: yaml.FullLoader
    :param tag_suffix: The name after the ``!func:`` tag
    :param type: str
    :param loader: A mapping node if function parameters given, scalar node if not
    :param type: Union[yaml.nodes.MappingNode, yaml.nodes.ScalarNode]
    """
    arg_map = None
    if isinstance(node, yaml.nodes.MappingNode):
        arg_map = (
            loader.construct_mapping(node, deep=True) if node.value else {}
        )
        node = yaml.nodes.ScalarNode(
            node.tag, "", node.start_mark, node.end_mark
        )
    func = loader.construct_python_name(tag_suffix, node)
    if arg_map:
        func = partial(func, **arg_map)
    return func

def get_loader() -> yaml.Loader:
    """Add constructors to PyYAML loader."""
    loader = yaml.FullLoader
    configs = [
        ClusterConfig.__name__,
    ]
    for config in configs:
        loader.add_constructor(f"!{config}", config_constructor)
    loader.add_multi_constructor("!func:", function_constructor)
    return loader

def parse_config(path: Path) -> BaseModel:
    """
    Parses a yaml configuration file into a config object.
    :param path: Path to configuration file
    :type path: pathlib.Path
    """

    with open(path, "rb") as f:
        yml_config = yaml.load(f, Loader=get_loader())
    return yml_config
